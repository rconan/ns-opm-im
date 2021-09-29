use crseo::{
    calibrations, ceo,
    dos::GmtOpticalSensorModel,
    from_opticals,
    shackhartmann::{Diffractive, Geometric, WavefrontSensorBuilder},
    Builder, Calibration, OpticalSensitivities, ShackHartmann, GMT, SH48, SOURCE,
};
use dosio::{ios, Dos};
use nalgebra::{DMatrix, DVector, SMatrix}; //, SVector
use osqp::{CscMatrix, Problem, Settings};
use serde::Deserialize;
use serde_pickle as pickle;
use skyangle::Conversion;
use std::time::Instant;
use std::{fs::File, io::BufReader};

// Shack-Hartmann WFS kind
type WfsType = Diffractive;
// Matrix type definitions
type MatrixNcxnc = SMatrix<f64, 271, 271>;
// Unable to use MatrixNsxnc data type (7360x271)
//type MatrixNsxnc = SMatrix<f64, 7360, 271>;
type DynMatrix = DMatrix<f64>;
//type VectorNs = SMatrix<f64, 7360, 1>;
type VectorNc = SMatrix<f64, 271, 1>;
//
const DEBUG_MSGS: bool = true;
// Number of bending modes (it can be retrieved from D)
const N_BM: u8 = 27;
// Ratio between cost J1 (WFS slope fitting) and J3 (control effort).
const J1_J3_RATIO: f64 = 10.0;
// Minimum value assigned to rho3 factor
const MIN_RHO3: f64 = 1.0e-6;

// AcO data structure
#[derive(Deserialize)]
struct QPData {
    #[serde(rename = "D")]
    dmat: Vec<f64>,
    #[serde(rename = "W2")]
    w2: Vec<f64>,
    #[serde(rename = "W3")]
    w3: Vec<f64>,
    #[serde(rename = "K")]
    k: f64,
    #[serde(rename = "wfsMask")]
    wfs_mask: Vec<Vec<bool>>,
    umin: Vec<f64>,
    umax: Vec<f64>,
    rm_mean_slopes: bool,
    #[serde(rename = "_Tu")]
    tu: Vec<f64>,
    rho_3: f64,
    end2end_ordering: bool,
}
#[derive(Deserialize)]
struct QP {
    #[serde(rename = "SHAcO_qp")]
    data: QPData,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // WFS Definition
    let n_aco_sensor = 3;
    let mut gosm_aco = GmtOpticalSensorModel::<ShackHartmann<WfsType>, SH48<WfsType>>::new(Some(
        SOURCE::new()
            .size(n_aco_sensor)
            .on_ring(6f32.from_arcmin())
            .fwhm(6f64),
    ))
    .gmt(GMT::new().m1("m1_eigen_modes", 329))
    .sensor(SH48::<WfsType>::new().n_sensor(n_aco_sensor))
    /*        .atmosphere(crseo::ATMOSPHERE::new().ray_tracing(
        26.,
        520,
        0.,
        25.,
        Some("ns-opm-im_atm.bin".to_string()),
        Some(8),
    ))*/
    .build()?;
    println!("valid lenslets: {}", &gosm_aco.sensor.n_valid_lenslet());

    let gmt_calib = ceo!(GMT, m1 = ["m1_eigen_modes", 27]); //GMT::new().m1("m1_eigen_modes", 27).build()?;
    let wfs_calib = SH48::<Geometric>::new().n_sensor(n_aco_sensor);
    let src_calib = wfs_calib
        .guide_stars(Some(
            SOURCE::new().size(n_aco_sensor).on_ring(6f32.from_arcmin()),
        ))
        .build()?;

    // WFS Calibration
    use calibrations::Mirror::*;
    use calibrations::Segment::*;
    let m1_rbm_sens: Vec<f32> = {
        println!("M1 RBM calibration");
        let mut calib = Calibration::new(&gmt_calib, &src_calib, wfs_calib.clone());
        let mut segments = vec![vec![Txyz(1e-6, None), Rxyz(1e-6, None)]; 6];
        segments.append(&mut vec![vec![Txyz(1e-6, None), Rxyz(1e-6, Some(0..2))]]);
        let now = Instant::now();
        calib.calibrate(
            vec![M1],
            segments,
            calibrations::ValidLensletCriteria::OtherSensor(&mut gosm_aco.sensor),
        );
        println!(
            "GMT WFS calibration [{}x{}] in {}s",
            &calib.n_data,
            &calib.n_mode,
            now.elapsed().as_secs()
        );
        calib.poke.into()
    };
    let m2_rbm_sens: Vec<f32> = {
        println!("M2 RBM calibration");
        let mut calib = Calibration::new(&gmt_calib, &src_calib, wfs_calib.clone());
        let mut segments = vec![vec![Txyz(1e-6, None), Rxyz(1e-6, None)]; 6];
        segments.append(&mut vec![vec![Txyz(1e-6, None), Rxyz(1e-6, Some(0..2))]]);
        let now = Instant::now();
        calib.calibrate(
            vec![M2],
            segments,
            calibrations::ValidLensletCriteria::OtherSensor(&mut gosm_aco.sensor),
        );
        println!(
            "GMT WFS calibration [{}x{}] in {}s",
            &calib.n_data,
            &calib.n_mode,
            now.elapsed().as_secs()
        );
        calib.poke.into()
    };
    let m1_bm_sens: Vec<f32> = {
        println!("M1 BM calibration");
        let mut calib = Calibration::new(&gmt_calib, &src_calib, wfs_calib.clone());
        let now = Instant::now();
        calib.calibrate(
            vec![M1MODES],
            vec![vec![Modes(1e-6, 0..27)]; 7],
            calibrations::ValidLensletCriteria::OtherSensor(&mut gosm_aco.sensor),
        );
        println!(
            "GMT WFS calibration [{}x{}] in {}s",
            &calib.n_data,
            &calib.n_mode,
            now.elapsed().as_secs()
        );
        calib.poke.into()
    };
    let dmat: Vec<_> = [m1_rbm_sens, m2_rbm_sens, m1_bm_sens]
        .iter()
        .cloned()
        .flatten()
        .map(|x| x as f64)
        .collect();

    // Import AcO data
    let qp: QP = {
        let file = File::open("SHAcO_qp_rhoP1e-3_kIp5.rs.pkl").unwrap();
        let rdr = BufReader::with_capacity(10_000, file);
        pickle::from_reader(rdr, Default::default()).unwrap()
    };
    println!("Number of lenslets:{}", qp.data.wfs_mask[0].len());

    // Handle loaded data
    let w2 = MatrixNcxnc::from_vec(qp.data.w2);
    let w3 = MatrixNcxnc::from_vec(qp.data.w3);

    let ns = dmat.len() / 271;
    //    println!("Valid lenslets:{}", qp.data.dmat.len() / 271);
    let d_wfs = DMatrix::from_column_slice(dmat.len() / 271, 271, &dmat);
    let d_wfs_svd = d_wfs.clone().svd(true, true);
    let sing_val_min = d_wfs_svd
        .singular_values
        .as_slice()
        .iter()
        .cloned()
        .fold(f64::INFINITY, |a, v| f64::min(a, v));
    let sing_val_max = d_wfs_svd
        .singular_values
        .as_slice()
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, |a, v| f64::max(a, v));
    println!("Condition number {}", sing_val_max / &sing_val_min);

    let d_t_w1_d = {
        let d_t_w1_d_dyn = d_wfs.tr_mul(&d_wfs);
        MatrixNcxnc::from_vec(d_t_w1_d_dyn.as_slice().to_vec())
    };

    // Extract the upper triangular elements of `P`
    let mut p_utri = {
        println!("rho_3:{}", qp.data.rho_3);
        let p = d_t_w1_d + w2 + w3.scale(qp.data.rho_3 * qp.data.k * qp.data.k);

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
        CscMatrix::from_column_iter_dense(p.nrows(), p.ncols(), p.as_slice().to_vec().into_iter())
            .into_upper_tri()
    };

    // Eigen modes and eigen coefs to forces matrix
    let file = File::open("../m1_eigen_modes.bin").unwrap();
    let data: Vec<(Vec<f64>, Vec<f64>)> = bincode::deserialize_from(file).unwrap();
    //    let (_, coefs2forces) = &data[sid - 1];
    let n_actuators = [335, 335, 335, 335, 335, 335, 306];
    let n_actuators_sum = n_actuators.iter().sum::<usize>();
    let fz = 10e-5;
    let coefs2forces: Vec<_> = data
        .iter()
        .zip(n_actuators.iter())
        .map(|((_, c2f), &n)| DMatrix::from_column_slice(n, c2f.len() / n, c2f) * fz)
        .inspect(|x| println!("{:?}", x.shape()))
        .collect();
    //    println!("coefs2forces: {:?}", coefs2force.shape());
    let tu_nrows = 41 * 2 + n_actuators_sum;
    let tu_ncols = 41 * 2 + 7 * N_BM as usize;
    let mut tu = DMatrix::<f64>::zeros(tu_nrows, tu_ncols);
    for (i, mut row) in tu.row_iter_mut().take(41).enumerate() {
        row[(i)] = 1f64;
    }
    for (i, mut row) in tu.row_iter_mut().skip(41).take(41).enumerate() {
        row[(i + 41)] = 1f64;
    }
    let mut n_skip_row = 82;
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

    // Remove S7Rz from T_u matrix
    // Indices to insert (or remove) S7Rz columns of matrix Tu
    /*
    let end2end_ordering = true;
    let i_m1_s7_rz: u8 = if end2end_ordering {
        41
    } else {
        ((12 + N_BM) * 6) + 5
    };
    let i_m2_s7_rz: u8 = if end2end_ordering {
        // Add 1 (+1) to delete
        82 + 1
    } else {
        ((12 + N_BM) * 6) + 10 + 1
    };
    let tu = DynMatrix::from_vec(273, 1228, qp.data.tu)
        .transpose()
        .remove_columns_at(&[i_m1_s7_rz.into(), i_m2_s7_rz.into()]);
    */
    // Inequality constraint matrix: lb <= a_in*u <= ub
    let a_in = {
        //println!("count nonzero: {}", qp.data.tu.iter().filter(|&n| *n != 0.0).count());
        let tus = tu.scale(qp.data.k);
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

    // ** Initialize u_ant vector --- mimic feedback
    //let u_ant = VectorNc::zeros();
    let u_ant = VectorNc::from_fn(|r, c| if r == c { 1.0e-6 } else { 1.0e-6 }) * 0f64;
    for i in 0..7 {
        println!("U_ant: {}", format!("{:.4e}", u_ant.get(i).unwrap()));
    }

    println!(
        "Number of valid lenslets:{}",
        gosm_aco.sensor.n_valid_lenslet()
    );
    let mut m1_modes_buf = vec![vec![0f64; 329]; 7];
    m1_modes_buf[0][0] = 1e-6;
    m1_modes_buf[3][5] = 1e-6;
    m1_modes_buf[6][3] = -1e-6;
    let m1_modes = ios!(M1modes(
        m1_modes_buf.into_iter().flatten().collect::<Vec<f64>>()
    ));
    let mut m1_rbm_buf = vec![vec![0f64; 6]; 7];
    m1_rbm_buf[1][1] = 1e-6;
    let m1_rbm = ios!(OSSM1Lcl(
        m1_rbm_buf.into_iter().flatten().collect::<Vec<f64>>()
    ));
    let mut m2_rbm_buf = vec![0f64; 42];
    m2_rbm_buf[3] = 1e-6;
    let m2_rbm = ios!(MCM2Lcl6D(m2_rbm_buf));
    let y_valid = gosm_aco
        .in_step_out(Some(vec![m1_rbm, m2_rbm, m1_modes]))?
        .and_then(|x| Into::<Option<Vec<f64>>>::into(x[0].clone()))
        .unwrap();
    println!("Number of valid slopes:{}", y_valid.len());
    println!("Valid slopes sum:{}", y_valid.iter().sum::<f64>());

    let y_vec = DVector::from_column_slice(&y_valid); //VectorNs::from_vec(y_valid);

    // QP linear term
    let mut q: Vec<f64> = (-y_vec.clone_owned().tr_mul(&d_wfs)
        - u_ant.tr_mul(&w3).scale(qp.data.rho_3 * qp.data.k))
    .as_slice()
    .to_vec();

    assert_eq!(271, q.len());

    // Update bounds to inequality constraints
    let umin = vec![f64::NEG_INFINITY; tu.nrows()];
    let umax = vec![f64::INFINITY; tu.nrows()];
    let tu_u_ant: Vec<f64> = (&tu * &u_ant).as_slice().to_vec();
    let lb: Vec<f64> = tu_u_ant
        .iter()
        .zip(umin.iter())
        .map(|(v, w)| w - v)
        .collect();
    let ub: Vec<f64> = tu_u_ant
        .iter()
        .zip(umax.iter())
        .map(|(v, w)| w - v)
        .collect();

    // QP settings
    let settings = Settings::default()
        .eps_abs(1.0e-8)
        .eps_rel(1.0e-6)
        .max_iter(500 * 271)
        .warm_start(true)
        .verbose(true);

    // Create an OSQP problem
    let mut prob =
        Problem::new(p_utri, &q, a_in, &lb, &ub, &settings).expect("Failed to setup AcO problem!");

    // Solve problem - 1st iteration
    let mut result = prob.solve();
    let mut c = result.x().expect("Failed to solve QP problem!");
    // Print the solution - Just first terms for verification
    for i in 0..7 {
        println!("{}", format!("{:.4e}", c[i]));
    }

    // Compute costs to set up the 2nd QP iteration
    let mut c_vec = VectorNc::from_vec(c.to_vec());
    let j_1na = {
        let epsilon = &y_vec - (&d_wfs * &c_vec);
        // Still need to account for W1
        epsilon.tr_mul(&epsilon)
    };
    // Control effort cost
    let j_3na = {
        let delta = c_vec.scale(qp.data.k) - u_ant;
        delta.tr_mul(&w3) * &delta
    };
    // nalgebra object to f64 scalar conversion
    let j_1 = j_1na.get(0).unwrap();
    let j_3 = j_3na.get(0).unwrap();

    if *j_3 != 0f64 {
        if DEBUG_MSGS {
            println!(
                "J1:{}J3:{}ratio:{}",
                format!("{:.4e} ", j_1),
                format!("{:.4e} ", j_3),
                format!("{:.4e} ", j_1 / (j_3 * qp.data.rho_3))
            );
        }

        let mut rho_3 = j_1 / (j_3 * J1_J3_RATIO);
        if rho_3 < MIN_RHO3 {
            rho_3 = MIN_RHO3
        };

        // Update QP P matrix
        p_utri = {
            println!("New rho_3:{}", format!("{:.4e}", rho_3));
            let p = d_t_w1_d + w2 + w3.scale(rho_3 * qp.data.k * qp.data.k);
            CscMatrix::from_column_iter_dense(
                p.nrows(),
                p.ncols(),
                p.as_slice().to_vec().into_iter(),
            )
            .into_upper_tri()
        };
        prob.update_P(p_utri);
        // Update QP linear term
        q = (-y_vec.clone_owned().tr_mul(&d_wfs) - u_ant.tr_mul(&w3).scale(rho_3 * qp.data.k))
            .as_slice()
            .to_vec();
        prob.update_lin_cost(&q);

        // Solve problem - 2nd iteration
        result = prob.solve();
        c = result.x().expect("Failed to solve QP problem!");
        c_vec = VectorNc::from_vec(c.to_vec());

        // Just for DEBUG
        if DEBUG_MSGS {
            // Print the solution - Just first terms for verification
            for i in 0..7 {
                println!("{}", format!("{:.4e}", c[i]));
            }

            let j_1na = {
                let epsilon = &y_vec - (&d_wfs * &c_vec);
                // Still need to account for W1
                epsilon.tr_mul(&epsilon)
            };
            // Control effort cost
            let j_3na = {
                let delta = c_vec.scale(qp.data.k) - u_ant;
                delta.tr_mul(&w3) * &delta
            };
            // nalgebra object to f64 scalar conversion
            let j_1 = j_1na.get(0).unwrap();
            let j_3 = j_3na.get(0).unwrap();

            println!(
                "J1:{}J3:{}ratio:{}",
                format!("{:.4e} ", j_1),
                format!("{:.4e} ", j_3),
                format!("{:.4e} ", j_1 / (j_3 * rho_3))
            );
        }
    }

    for k in 0..7 {
        println!(
            "M1 S{} RBM Solution: {:#?}",
            k + 1,
            c_vec
                .rows(0 + k * 6, if k == 6 { 5 } else { 6 })
                .as_slice()
                .iter()
                .map(|x| x * 1e6)
                .collect::<Vec<f64>>()
        );
        println!(
            "M2 S{} RBM Solution: {:#?}",
            k + 1,
            c_vec
                .rows(41 + k * 6, if k == 6 { 5 } else { 6 })
                .as_slice()
                .iter()
                .map(|x| x * 1e6)
                .collect::<Vec<f64>>()
        );
        println!(
            "M1 S{} BM Solution: {:#?}",
            k + 1,
            c_vec
                .rows(82 + k * 27, 27)
                .as_slice()
                .iter()
                .map(|x| x * 1e6)
                .collect::<Vec<f64>>()
        );
    }
    Ok(())
}
