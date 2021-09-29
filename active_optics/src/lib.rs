use dosio::{ios, DOSIOSError, Dos, IOVec, IO};
use nalgebra as na;
use nalgebra::{DMatrix, DVector, SMatrix}; //, SVector
use osqp::{CscMatrix, Problem, Settings};
use serde::Deserialize;
use serde_pickle as pickle;
use std::{fs::File, io::BufReader};

// Matrix type definitions
type DynMatrix = DMatrix<f64>;
// Ratio between cost J1 (WFS slope fitting) and J3 (control effort).
const J1_J3_RATIO: f64 = 10.0;
// Minimum value assigned to rho3 factor
const MIN_RHO3: f64 = 1.0e-6;

// AcO data structure
#[derive(Deserialize)]
struct QPData {
    #[serde(rename = "W2")]
    w2: Vec<f64>,
    #[serde(rename = "W3")]
    w3: Vec<f64>,
    #[serde(rename = "K")]
    k: f64,
    rho_3: f64,
}
#[derive(Deserialize)]
pub struct QP<'a, const M1_RBM: usize, const M2_RBM: usize, const M1_BM: usize, const N_MODE: usize>
{
    #[serde(rename = "SHAcO_qp")]
    data: QPData,
    /// calibration matrix (column wise) as (data,n_cols)
    #[serde(skip)]
    dmat: (&'a [f64], usize),
    /// segment bending modes coefficients to segment actuators forces  (column wise) as ([data],[n_rows])
    #[serde(skip)]
    coefs2forces: (&'a [Vec<f64>], Vec<usize>),
}

impl<'a, const M1_RBM: usize, const M2_RBM: usize, const M1_BM: usize, const N_MODE: usize>
    QP<'a, M1_RBM, M2_RBM, M1_BM, N_MODE>
{
    pub fn new(
        qp_filename: &str,
        calib_matrix: (&'a [f64], usize),
        coefs2forces: (&'a [Vec<f64>], Vec<usize>),
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(qp_filename)?;
        let rdr = BufReader::with_capacity(10_000, file);
        let mut this: Self = pickle::from_reader(rdr, Default::default())?;
        this.dmat = calib_matrix;
        this.coefs2forces = coefs2forces;
        Ok(this)
    }
    fn modal2nodal(&self) -> DynMatrix {
        let (data, n_actuators) = &self.coefs2forces;
        let n_actuators_sum = n_actuators.iter().sum::<usize>();
        let fz = 10e-5;
        let coefs2forces: Vec<_> = data
            .iter()
            .zip(n_actuators.iter())
            .map(|(c2f, &n)| DMatrix::from_column_slice(n, c2f.len() / n, c2f) * fz)
            .inspect(|x| println!("{:?}", x.shape()))
            .collect();
        //    println!("coefs2forces: {:?}", coefs2force.shape());
        let tu_nrows = M1_RBM * M2_RBM + n_actuators_sum;
        let tu_ncols = N_MODE;
        let mut tu = DMatrix::<f64>::zeros(tu_nrows, tu_ncols);
        for (i, mut row) in tu.row_iter_mut().take(M1_RBM).enumerate() {
            row[(i)] = 1f64;
        }
        for (i, mut row) in tu.row_iter_mut().skip(M1_RBM).take(M2_RBM).enumerate() {
            row[(i + M1_RBM)] = 1f64;
        }
        let mut n_skip_row = M1_RBM + M2_RBM;
        for c2f in coefs2forces {
            for (mut tu_row, c2f_row) in tu.row_iter_mut().skip(n_skip_row).zip(c2f.row_iter()) {
                tu_row
                    .iter_mut()
                    .zip(c2f_row.iter())
                    .for_each(|(tu, &c2f)| {
                        *tu = c2f;
                    });
                n_skip_row += c2f.nrows();
            }
        }
        tu
    }
    pub fn build(self) -> ActiveOptics<N_MODE> {
        let (dmat, n_mode) = self.dmat;
        assert!(n_mode == N_MODE);
        // W2 and W3
        let w2 = SMatrix::<f64, N_MODE, N_MODE>::from_column_slice(&self.data.w2);
        let w3 = SMatrix::<f64, N_MODE, N_MODE>::from_column_slice(&self.data.w3);
        // W1
        let d_wfs = DMatrix::from_column_slice(dmat.len() / n_mode, n_mode, &dmat);
        let d_t_w1_d = {
            let d_t_w1_d_dyn = d_wfs.tr_mul(&d_wfs);
            SMatrix::<f64, N_MODE, N_MODE>::from_vec(d_t_w1_d_dyn.as_slice().to_vec())
        };
        // Extract the upper triangular elements of `P`
        let p_utri = {
            println!("rho_3:{}", self.data.rho_3);
            let p = d_t_w1_d + w2 + w3.scale(self.data.rho_3 * self.data.k * self.data.k);

            // Check matrix density
            let p_nnz = p
                .as_slice()
                .iter()
                .filter_map(|&p| if p != 0.0 { Some(1.0) } else { None })
                .sum::<f64>();
            println!(
                "P: {:?}, density: {}%",
                p.shape(),
                100. * p_nnz / (p.ncols() * p.nrows()) as f64
            );
            CscMatrix::from_column_iter_dense(
                p.nrows(),
                p.ncols(),
                p.as_slice().to_vec().into_iter(),
            )
            .into_upper_tri()
        };

        let tu = self.modal2nodal();

        let a_in = {
            //println!("count nonzero: {}", self.data.tu.iter().filter(|&n| *n != 0.0).count());
            let tus = tu.scale(self.data.k);
            let tu_nnz = tus.as_slice().iter().fold(0.0, |mut s, p| {
                if *p != 0.0 {
                    s += 1.0;
                };
                s
            });
            println!(
                "Tu: {:?}, nnz: {}, density: {:.0}%",
                tus.shape(),
                tu_nnz,
                100. * tu_nnz / (tus.ncols() * tus.nrows()) as f64
            );

            println!("Number of Tu cols:{}", tu.ncols());

            CscMatrix::from(
                &tus.row_iter()
                    .map(|x| x.clone_owned().as_slice().to_vec())
                    .collect::<Vec<Vec<f64>>>(),
            )
        };
        println!("a_in: {:?}", (a_in.nrows, a_in.ncols));

        // QP linear term
        let q: Vec<f64> = vec![0f64; N_MODE];

        // Inequality constraints
        let umin = vec![f64::NEG_INFINITY; tu.nrows()];
        let umax = vec![f64::INFINITY; tu.nrows()];

        // QP settings
        let settings = Settings::default()
            .eps_abs(1.0e-8)
            .eps_rel(1.0e-6)
            .max_iter(500 * 271)
            .warm_start(true)
            .verbose(true);

        // Create an OSQP problem
        let prob = Problem::new(p_utri, &q, a_in, &umin, &umax, &settings)
            .expect("Failed to setup AcO problem!");

        ActiveOptics {
            prob,
            u: Vec::with_capacity(84 + 7 * M1_BM),
            y_valid: Vec::with_capacity(d_wfs.nrows()),
            d_wfs,
            u_ant: SMatrix::zeros(),
            d_t_w1_d,
            w2,
            w3,
            rho_3: self.data.rho_3,
            k: self.data.k,
            umin: vec![f64::NEG_INFINITY; tu.nrows()],
            umax: vec![f64::INFINITY; tu.nrows()],
            tu,
        }
    }
}

pub struct ActiveOptics<const N_MODE: usize> {
    prob: Problem,
    d_wfs: DynMatrix,
    u_ant: SMatrix<f64, N_MODE, 1>,
    u: Vec<f64>,
    y_valid: Vec<f64>,
    d_t_w1_d: na::Matrix<
        f64,
        na::Const<N_MODE>,
        na::Const<N_MODE>,
        na::ArrayStorage<f64, N_MODE, N_MODE>,
    >,
    w2: SMatrix<f64, N_MODE, N_MODE>,
    w3: SMatrix<f64, N_MODE, N_MODE>,
    rho_3: f64,
    k: f64,
    umin: Vec<f64>,
    umax: Vec<f64>,
    tu: DynMatrix,
}

impl<const N_MODE: usize> Iterator for ActiveOptics<N_MODE> {
    type Item = ();

    fn next(&mut self) -> Option<Self::Item> {
        let y_vec = DVector::from_column_slice(&self.y_valid); //VectorNs::from_vec(y_valid);
                                                               // QP linear term                                                               // QP linear term
        let mut q: Vec<f64> = (-y_vec.clone_owned().tr_mul(&self.d_wfs)
            - self.u_ant.tr_mul(&self.w3).scale(self.rho_3 * self.k))
        .as_slice()
        .to_vec();
        self.prob.update_lin_cost(&q);
        // Update bounds to inequality constraints
        let tu_u_ant: Vec<f64> = (&self.tu * &self.u_ant).as_slice().to_vec();
        let lb: Vec<f64> = tu_u_ant
            .iter()
            .zip(self.umin.iter())
            .map(|(v, w)| w - v)
            .collect();
        let ub: Vec<f64> = tu_u_ant
            .iter()
            .zip(self.umax.iter())
            .map(|(v, w)| w - v)
            .collect();
        self.prob.update_bounds(&lb, &ub);
        let mut result = self.prob.solve();
        let mut c = result.x().expect("Failed to solve QP problem!");
        // Compute costs to set up the 2nd QP iteration
        let c_vec = SMatrix::<f64, N_MODE, 1>::from_vec(c.to_vec());
        let j_1na = {
            let epsilon = &y_vec - (&self.d_wfs * &c_vec);
            // Still need to account for W1
            epsilon.tr_mul(&epsilon)
        };
        // Control effort cost
        let j_3na = {
            let delta = c_vec.scale(self.k) - self.u_ant;
            delta.tr_mul(&self.w3) * &delta
        };
        // nalgebra object to f64 scalar conversion
        let j_1 = j_1na.get(0).unwrap();
        let j_3 = j_3na.get(0).unwrap();
        if *j_3 != 0f64 {
            let mut rho_3 = j_1 / (j_3 * J1_J3_RATIO);
            if rho_3 < MIN_RHO3 {
                rho_3 = MIN_RHO3
            };

            // Update QP P matrix
            let p_utri = {
                //                println!("New rho_3:{}", format!("{:.4e}", rho_3));
                let p = self.d_t_w1_d + self.w2 + self.w3.scale(rho_3 * self.k * self.k);
                CscMatrix::from_column_iter_dense(
                    p.nrows(),
                    p.ncols(),
                    p.as_slice().to_vec().into_iter(),
                )
                .into_upper_tri()
            };
            self.prob.update_P(p_utri);
            // Update QP linear term
            q = (-y_vec.clone_owned().tr_mul(&self.d_wfs)
                - self.u_ant.tr_mul(&self.w3).scale(rho_3 * self.k))
            .as_slice()
            .to_vec();
            self.prob.update_lin_cost(&q);

            // Solve problem - 2nd iteration
            result = self.prob.solve();
            c = result.x().expect("Failed to solve QP problem!");
            //            c_vec = SMatrix::<f64, N_MODE, 1>::from_vec(c.to_vec());
            let k = self.k;
            self.u[..41]
                .iter_mut()
                .zip(&c[..41])
                .for_each(|(u, c)| *u -= k * c);
            self.u[42..83]
                .iter_mut()
                .zip(&c[41..82])
                .for_each(|(u, c)| *u -= k * c);
            self.u[84..]
                .iter_mut()
                .zip(&c[82..])
                .for_each(|(u, c)| *u -= k * c);
        }
        Some(())
    }
}

impl<const N_MODE: usize> Dos for ActiveOptics<N_MODE> {
    type Input = Vec<f64>;
    type Output = Vec<f64>;

    fn outputs(&mut self) -> Option<Vec<IO<Self::Output>>> {
        let segment_bm: Vec<_> = self.u[84..].chunks(27).collect();
        Some(ios!(
            M1RBMcmd(self.u[..42].to_vec()),
            M2poscmd(self.u[42..84].to_vec()),
            M1S1BMcmd(segment_bm[0].to_vec()),
            M1S2BMcmd(segment_bm[1].to_vec()),
            M1S3BMcmd(segment_bm[2].to_vec()),
            M1S4BMcmd(segment_bm[3].to_vec()),
            M1S5BMcmd(segment_bm[4].to_vec()),
            M1S6BMcmd(segment_bm[5].to_vec()),
            M1S7BMcmd(segment_bm[6].to_vec())
        ))
    }

    fn inputs(
        &mut self,
        data: Option<Vec<IO<Self::Input>>>,
    ) -> Result<&mut Self, dosio::DOSIOSError> {
        match data {
            Some(mut data) => {
                if let Some(IO::SensorData { data: Some(value) }) = data.pop_this(ios!(SensorData))
                {
                    self.y_valid = value;
                    Ok(self)
                } else {
                    Err(DOSIOSError::Inputs(
                        "ActiveOptics SensorData not found".into(),
                    ))
                }
            }
            None => Err(DOSIOSError::Inputs(
                "None data passed to Active Optics".into(),
            )),
        }
    }
}
