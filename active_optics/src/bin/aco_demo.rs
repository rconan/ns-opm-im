//
// -----------------------------------------------------------------
// First version of an AcO standalone implementation
// -----------------------------------------------------------------
// It does not implement the mean slop removal feature.
// August, 2021
//
// R. Romano & R. Conan
//
// .............................................................................

use nalgebra::{DMatrix, SMatrix}; //, SVector
                                  //use nalgebra_sparse::csc::CscMatrix as naCSC;
use osqp::{CscMatrix, Problem, Settings};
use serde::Deserialize;
use serde_pickle as pickle;
use std::{fs::File, io::BufReader};
// .............................................................................

// Matrix type definitions
type MatrixNcxnc = SMatrix<f64, 271, 271>;
// Unable to use MatrixNsxnc data type (7360x271)
//type MatrixNsxnc = SMatrix<f64, 7360, 271>;
type DynMatrix = DMatrix<f64>;
type VectorNs = SMatrix<f64, 7360, 1>;
type VectorNc = SMatrix<f64, 271, 1>;
// .............................................................................

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

// wfs48x48 sample structure
#[derive(Deserialize)]
struct WFSData {
    wfsdata: Vec<f32>,
}
// .............................................................................

//
const DEBUG_MSGS: bool = true;
// Number of bending modes (it can be retrieved from D)
const N_BM: u8 = 27;
// Ratio between cost J1 (WFS slope fitting) and J3 (control effort).
const J1_J3_RATIO: f64 = 10.0;
// Minimum value assigned to rho3 factor
const MIN_RHO3: f64 = 1.0e-6;

//
fn get_valid_y(s_struct: WFSData, wfs_mask: Vec<Vec<bool>>) -> Vec<f64> {
    let mut y_valid = Vec::new();
    for seg_mask in wfs_mask.iter() {
        // Commands to take indices of true elements:
        // https://codereview.stackexchange.com/questions/159652/indices-of-true-values
        let valid_l_idxs: Vec<_> = seg_mask
            .iter()
            .enumerate()
            .filter(|&(_, &value)| value)
            .map(|(index, _)| index)
            .collect();
        //println!("...{}", valid_l_idxs.len());

        for iv in valid_l_idxs {
            y_valid.push(s_struct.wfsdata[iv] as f64);
        }
        // ::retain may be an alternative
        // https://doc.rust-lang.org/std/vec/struct.Vec.html#method.retain
    }
    return y_valid;
}

// .............................................................................

fn main() {
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
    let ns = qp.data.dmat.len() / 271;
    println!("Valid lenslets:{}", qp.data.dmat.len() / 271);
    let d_wfs = DynMatrix::from_vec(271, 7360, qp.data.dmat).transpose();

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

    // Remove S7Rz from T_u matrix
    // Indices to insert (or remove) S7Rz columns of matrix Tu
    let i_m1_s7_rz: u8 = if qp.data.end2end_ordering {
        41
    } else {
        ((12 + N_BM) * 6) + 5
    };
    let i_m2_s7_rz: u8 = if qp.data.end2end_ordering {
        // Add 1 (+1) to delete
        82 + 1
    } else {
        ((12 + N_BM) * 6) + 10 + 1
    };
    let tu = DynMatrix::from_vec(273, 1228, qp.data.tu)
        .transpose()
        .remove_columns_at(&[i_m1_s7_rz.into(), i_m2_s7_rz.into()]);

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

    // ** Initialize u_ant vector --- mimic feedback
    //let u_ant = VectorNc::zeros();
    let u_ant = VectorNc::from_fn(|r, c| if r == c { 1.0e-6 } else { 1.0e-6 });
    for i in 0..7 {
        println!("{}", format!("{:.4e}", u_ant.get(i).unwrap()));
    }

    // Import AcO data --- mimic feedback
    let s_struct: WFSData = {
        let file = File::open("wfs48x48sample.rs.pkl").unwrap();
        let rdr = BufReader::with_capacity(100_000, file);
        pickle::from_reader(rdr, Default::default()).unwrap()
    };
    println!("Number of lenslets:{}", s_struct.wfsdata.len());
    let y_valid = get_valid_y(s_struct, qp.data.wfs_mask);
    assert_eq!(ns, y_valid.len()); // WFS meas dimension check
                                   // Save as ns-dimensional nalgebra vector
    let y_vec = VectorNs::from_vec(y_valid);

    // QP linear term
    let mut q: Vec<f64> = (-y_vec.clone_owned().tr_mul(&d_wfs)
        - u_ant.tr_mul(&w3).scale(qp.data.rho_3 * qp.data.k))
    .as_slice()
    .to_vec();
    assert_eq!(271, q.len());

    // Update bounds to inequality constraints
    let tu_u_ant: Vec<f64> = (&tu * &u_ant).as_slice().to_vec();
    let lb: Vec<f64> = tu_u_ant
        .iter()
        .zip(qp.data.umin.iter())
        .map(|(v, w)| w - v)
        .collect();
    let ub: Vec<f64> = tu_u_ant
        .iter()
        .zip(qp.data.umax.iter())
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
        CscMatrix::from_column_iter_dense(p.nrows(), p.ncols(), p.as_slice().to_vec().into_iter())
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
