//! Active Optics Control Algorithm
//!
//! Implementation of the GMT Active Optics Control Algorithm
//! as described in "Natural Seeing Control Algorithms Description
//! Document" (GMT-DOC-04426)

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

/// Data structure for the quadratic programming algorithm
#[derive(Deserialize)]
struct QPData {
    /// Controllable mode regularization matrix
    #[serde(rename = "W2")]
    w2: Vec<f64>,
    /// Control balance weighting matrix
    #[serde(rename = "W3")]
    w3: Vec<f64>,
    /// Controller gain
    #[serde(rename = "K")]
    k: f64,
    /// Objective function factor
    rho_3: f64,
}
/// Quadratic programming stucture
///
/// It requires 4 generic constants:
///  - `M1_RBM`: the number of controlled M1 segment rigid body motions (at most 41 as S7 Rz is a blind mode)
///  - `M2_RBM`: the number of controlled M2 segment rigid body motions (at most 41 as S7 Rz is a blind mode)
///  - `M1_BM` : the number of controlled M1 segment eigen modes
///  - `N_MODE` = M1_RBM + M2_RBM + 7 * M1_BM
#[derive(Deserialize)]
pub struct QP<'a, const M1_RBM: usize, const M2_RBM: usize, const M1_BM: usize, const N_MODE: usize>
{
    /// Quadratic programming data
    #[serde(rename = "SHAcO_qp")]
    data: QPData,
    /// calibration matrix (column wise) as (data,n_cols)
    #[serde(skip)]
    dmat: (&'a [f64], usize),
    /// segment bending modes coefficients to segment actuators forces  (column wise) as ([data],[n_rows])
    #[serde(skip)]
    coefs2forces: (&'a [Vec<f64>], Vec<usize>),
    /// OSQP verbosity
    #[serde(skip)]
    verbose: bool,
    #[serde(skip)]
    /// convert bending modes coefficients to forces if true
    m1_actuator_forces_outputs: bool,
    /// OSQP convergence tolerances (absolute=1e-8,relative=1e-6)
    #[serde(skip)]
    convergence_tolerances: (f64, f64),
}

impl<'a, const M1_RBM: usize, const M2_RBM: usize, const M1_BM: usize, const N_MODE: usize>
    QP<'a, M1_RBM, M2_RBM, M1_BM, N_MODE>
{
    /// Creates a new quadratic programming object
    pub fn new(
        qp_filename: &str,
        calib_matrix: (&'a [f64], usize),
        coefs2forces: (&'a [Vec<f64>], Vec<usize>),
    ) -> Result<Self, Box<dyn std::error::Error>> {
        assert!(
            M1_RBM + M2_RBM + 7 * M1_BM == N_MODE,
            "The number of modes {} do not match the expected value {} (M1_RBM + M2_RBM + 7 * M1_BM)x",
            N_MODE, M1_RBM + M2_RBM + 7 * M1_BM
        );
        let file = File::open(qp_filename)?;
        let rdr = BufReader::with_capacity(10_000, file);
        let this: Self = pickle::from_reader(rdr, Default::default())?;
        Ok(Self {
            dmat: calib_matrix,
            coefs2forces,
            verbose: false,
            m1_actuator_forces_outputs: false,
            convergence_tolerances: (1.0e-8, 1.0e-6),
            ..this
        })
    }
    /// Sets OSQP verbosity
    pub fn verbose(self) -> Self {
        Self {
            verbose: true,
            ..self
        }
    }
    /// Outputs the forces of M1 actuators instead of M1 bending modes coefficients
    pub fn as_m1_actuator_forces(self) -> Self {
        Self {
            m1_actuator_forces_outputs: true,
            ..self
        }
    }
    /// Computes the control balance weighting matrix
    fn balance_weighting(&self) -> DynMatrix {
        let (data, n_actuators) = &self.coefs2forces;
        let n_actuators_sum = n_actuators.iter().sum::<usize>();
        let fz = 10e-5;
        let coefs2forces: Vec<_> = data
            .iter()
            .zip(n_actuators.iter())
            .map(|(c2f, &n)| DMatrix::from_column_slice(n, c2f.len() / n, c2f) * fz)
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
    /// Sets OSQP convergence tolerances: (absolute,relative)
    pub fn convergence_tolerances(self, convergence_tolerances: (f64, f64)) -> Self {
        Self {
            convergence_tolerances,
            ..self
        }
    }
    /// Builds the quadratic programming problem
    pub fn build(self) -> ActiveOptics<M1_RBM, M2_RBM, M1_BM, N_MODE> {
        let (dmat, n_mode) = self.dmat;
        assert!(n_mode == N_MODE,"The number of columns ({}) of the calibration matrix do not match the number of modes ({})",n_mode,N_MODE);
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
            let p = d_t_w1_d + w2 + w3.scale(self.data.rho_3 * self.data.k * self.data.k);
            CscMatrix::from_column_iter_dense(
                p.nrows(),
                p.ncols(),
                p.as_slice().to_vec().into_iter(),
            )
            .into_upper_tri()
        };

        let tu = self.balance_weighting();

        let a_in = {
            let tus = tu.scale(self.data.k);
            CscMatrix::from(
                &tus.row_iter()
                    .map(|x| x.clone_owned().as_slice().to_vec())
                    .collect::<Vec<Vec<f64>>>(),
            )
        };

        // QP linear term
        let q: Vec<f64> = vec![0f64; N_MODE];

        // Inequality constraints
        let umin = vec![f64::NEG_INFINITY; tu.nrows()];
        let umax = vec![f64::INFINITY; tu.nrows()];

        // QP settings
        let settings = Settings::default()
            .eps_abs(self.convergence_tolerances.0)
            .eps_rel(self.convergence_tolerances.1)
            .max_iter(500 * 271)
            .warm_start(true)
            .verbose(self.verbose);

        // Create an OSQP problem
        let prob = Problem::new(p_utri, &q, a_in, &umin, &umax, &settings)
            .expect("Failed to setup AcO problem!");

        ActiveOptics {
            prob,
            u: vec![0f64; 84 + 7 * M1_BM],
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
            coefs2forces: self.m1_actuator_forces_outputs.then(|| {
                let (data, n_actuators) = &self.coefs2forces;
                data.iter()
                    .zip(n_actuators)
                    .map(|(c2f, &n)| {
                        DMatrix::from_column_slice(n, c2f.len() / n, c2f)
                            .columns(0, M1_BM)
                            .into_owned()
                    })
                    .collect()
            }),
        }
    }
}

pub struct ActiveOptics<
    const M1_RBM: usize,
    const M2_RBM: usize,
    const M1_BM: usize,
    const N_MODE: usize,
> {
    /// Quadratic programming problem
    prob: Problem,
    /// Calibration matrix
    d_wfs: DynMatrix,
    /// Previous quadratic programming solution
    u_ant: SMatrix<f64, N_MODE, 1>,
    /// Current quadratic programming solution
    u: Vec<f64>,
    /// Wavefront sensor data
    y_valid: Vec<f64>,
    /// Wavefront error weighting matrix
    d_t_w1_d: na::Matrix<
        f64,
        na::Const<N_MODE>,
        na::Const<N_MODE>,
        na::ArrayStorage<f64, N_MODE, N_MODE>,
    >,
    /// Controllable mode regularization matrix    
    w2: SMatrix<f64, N_MODE, N_MODE>,
    /// Control balance weighting matrix
    w3: SMatrix<f64, N_MODE, N_MODE>,
    /// Objective function factor
    rho_3: f64,
    /// Controller gain
    k: f64,
    /// QP solution lower bound
    umin: Vec<f64>,
    /// QP solution upper bound
    umax: Vec<f64>,
    /// Control balance weighting matrix
    tu: DynMatrix,
    /// segment bending modes coefficients to segment actuators forces  (column wise) as ([data],[n_rows])
    coefs2forces: Option<Vec<DynMatrix>>,
}
impl<const M1_RBM: usize, const M2_RBM: usize, const M1_BM: usize, const N_MODE: usize>
    ActiveOptics<M1_RBM, M2_RBM, M1_BM, N_MODE>
{
    pub fn controller_gain(&self) -> f64 {
        self.k
    }
}

impl<const M1_RBM: usize, const M2_RBM: usize, const M1_BM: usize, const N_MODE: usize> Iterator
    for ActiveOptics<M1_RBM, M2_RBM, M1_BM, N_MODE>
{
    type Item = ();
    /// Updates the quadratic programming problem and computes the new solution
    fn next(&mut self) -> Option<Self::Item> {
        let y_vec = DVector::from_column_slice(&self.y_valid); //VectorNs::from_vec(y_valid);

        self.u_ant
            .iter_mut()
            .zip(&self.u)
            .take(M1_RBM)
            .for_each(|(u, &v)| *u = v);
        self.u_ant
            .iter_mut()
            .skip(M1_RBM)
            .zip(&self.u[42..])
            .take(M2_RBM)
            .for_each(|(u, &v)| *u = v);
        self.u_ant
            .iter_mut()
            .skip(M1_RBM + M2_RBM)
            .zip(&self.u[84..])
            .for_each(|(u, &v)| *u = v);

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
        let mut c = match result.x() {
            Some(x) => x,
            None => {
                pickle::to_writer(
                    &mut File::create("OSQP_log.pkl").unwrap(),
                    &(self.y_valid.clone(), self.u_ant.as_slice().to_vec()),
                    Default::default(),
                )
                .unwrap();
                panic!("Failed to solve QP problem!");
            }
        }; //.expect("Failed to solve QP problem!");
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
        //println!(" ===>>> J3: {:}:, J1: {:},RHO_3: {:}", j_3, j_1, self.rho_3);
        if *j_3 != 0f64 {
            self.rho_3 = j_1 / (j_3 * J1_J3_RATIO);
            if self.rho_3 < MIN_RHO3 {
                self.rho_3 = MIN_RHO3
            };

            // Update QP P matrix
            let p_utri = {
                //println!("New rho_3:{}", format!("{:.4e}", self.rho_3));
                let p = self.d_t_w1_d + self.w2 + self.w3.scale(self.rho_3 * self.k * self.k);
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
                - self.u_ant.tr_mul(&self.w3).scale(self.rho_3 * self.k))
            .as_slice()
            .to_vec();
            self.prob.update_lin_cost(&q);

            // Solve problem - 2nd iteration
            result = self.prob.solve();
            c = match result.x() {
                Some(x) => x,
                None => {
                    pickle::to_writer(
                        &mut File::create("OSQP_log.pkl").unwrap(),
                        &(self.y_valid.clone(), self.u_ant.as_slice().to_vec()),
                        Default::default(),
                    )
                    .unwrap();
                    panic!("Failed to solve QP problem!");
                }
            }; //.expect("Failed to solve QP problem!");
        }
        // Control action
        let k = self.k;
        self.u
            .iter_mut()
            .zip(&c[..M1_RBM])
            .for_each(|(u, c)| *u -= k * c); // u = u - k * c
        self.u[42..]
            .iter_mut()
            .zip(&c[M1_RBM..M1_RBM + M2_RBM])
            .for_each(|(u, c)| *u -= k * c);
        self.u[84..]
            .iter_mut()
            .zip(&c[M1_RBM + M2_RBM..])
            .for_each(|(u, c)| *u -= k * c);
        Some(())
    }
}

impl<const M1_RBM: usize, const M2_RBM: usize, const M1_BM: usize, const N_MODE: usize> Dos
    for ActiveOptics<M1_RBM, M2_RBM, M1_BM, N_MODE>
{
    type Input = Vec<f64>;
    type Output = Vec<f64>;

    fn outputs(&mut self) -> Option<Vec<IO<Self::Output>>> {
        match &self.coefs2forces {
            Some(c2f) => {
                let mut segment_bm = self.u[84..]
                    .chunks(M1_BM)
                    .map(|u| na::DVector::from_column_slice(u))
                    .zip(c2f)
                    .map(|(u, c2f)| c2f * u);
                Some(ios!(
                    M1RBMcmd(self.u[..42].to_vec()),
                    M2poscmd(self.u[42..84].to_vec()),
                    M1S1BMcmd(segment_bm.next().unwrap().as_slice().to_vec()),
                    M1S2BMcmd(segment_bm.next().unwrap().as_slice().to_vec()),
                    M1S3BMcmd(segment_bm.next().unwrap().as_slice().to_vec()),
                    M1S4BMcmd(segment_bm.next().unwrap().as_slice().to_vec()),
                    M1S5BMcmd(segment_bm.next().unwrap().as_slice().to_vec()),
                    M1S6BMcmd(segment_bm.next().unwrap().as_slice().to_vec()),
                    M1S7BMcmd(segment_bm.next().unwrap().as_slice().to_vec())
                ))
            }
            None => {
                let mut segment_bm = self.u[84..].chunks(M1_BM);
                Some(ios!(
                    M1RBMcmd(self.u[..42].to_vec()),
                    M2poscmd(self.u[42..84].to_vec()),
                    M1S1BMcmd(segment_bm.next().unwrap().to_vec()),
                    M1S2BMcmd(segment_bm.next().unwrap().to_vec()),
                    M1S3BMcmd(segment_bm.next().unwrap().to_vec()),
                    M1S4BMcmd(segment_bm.next().unwrap().to_vec()),
                    M1S5BMcmd(segment_bm.next().unwrap().to_vec()),
                    M1S6BMcmd(segment_bm.next().unwrap().to_vec()),
                    M1S7BMcmd(segment_bm.next().unwrap().to_vec())
                ))
            }
        }
    }

    fn inputs(
        &mut self,
        data: Option<Vec<IO<Self::Input>>>,
    ) -> Result<&mut Self, dosio::DOSIOSError> {
        match data {
            Some(mut data) => match data.pop_this(ios!(SensorData)) {
                Some(IO::SensorData { data: Some(value) }) if value.len() == self.d_wfs.nrows() => {
                    self.y_valid = value;
                    Ok(self)
                }
                _ => Err(DOSIOSError::Inputs(
                    "No suitable ActiveOptics SensorData found, either the data is None or its size do not match the calibration matrix rows".into(),
                )),
            },
            None => Err(DOSIOSError::Inputs(
                "None data passed to Active Optics".into(),
            )),
        }
    }
}
