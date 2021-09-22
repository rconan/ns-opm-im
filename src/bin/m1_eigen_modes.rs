use complot::{tri, Config};
use fem::FEM;
use geotrans::{Quaternion, Vector};
use indicatif::ProgressBar;
use nalgebra as na;
use rand::prelude::*;
use rayon::prelude::*;
use serde_pickle as pkl;
use spade::{delaunay::FloatDelaunayTriangulation, HasPosition};
use std::error::Error;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use triangle_rs as triangle;

type Matrix =
    na::Matrix<f64, na::Dynamic, na::Dynamic, na::VecStorage<f64, na::Dynamic, na::Dynamic>>;
type U_BM2F = (Matrix, Matrix, Vec<f64>);
fn heatmap(filename: &str, map: &[f64], nodes: &[f64], scaling: Option<f64>) {
    let delaunay = triangle::Builder::new()
        .set_tri_points(
            nodes
                .clone()
                .chunks(3)
                .map(|x| x[..2].to_vec())
                .flatten()
                .collect::<Vec<f64>>(),
        )
        .set_switches("Q")
        .build();

    let cells: Vec<f64> = delaunay
        .triangle_iter()
        .map(|t| t.iter().fold(0., |a, &i| a + map[i] / 3.))
        .collect();

    let data = delaunay
        .triangle_vertex_iter()
        .zip(cells.into_iter().map(|x| x * scaling.unwrap_or(1f64)));
    tri::Heatmap::from((data, Some(Config::new().filename(filename))));
}

fn m1_actuators_2_rbm(sid: usize) -> Result<Matrix, Box<dyn Error>> {
    let mut fem =
        FEM::from_pickle("data/20210802_0755_MT_mount_v202104_FSM/static_reduction_model.73.pkl")?;
    let n_io = (fem.n_inputs(), fem.n_outputs());
    //    println!("{}", fem);
    fem.keep_inputs(&[sid - 1]);
    fem.keep_outputs_by(&[25], |x| x.descriptions.contains(&format!("M1-S{}", sid)));

    let gain = fem.reduced_static_gain(n_io).unwrap();
    println!("M1 S{} actuactors / RBM gain: [{:?}]", sid, &gain.shape());
    let gain_svd = gain.svd(true, true);
    println!("Singular values:");
    println!("{:#?}", &gain_svd.singular_values.as_slice());
    Ok(gain_svd.v_t.unwrap())
}

fn apply_force(
    sid: usize,
    force: na::Matrix<
        f64,
        na::Dynamic,
        na::Const<1_usize>,
        na::VecStorage<f64, na::Dynamic, na::Const<1_usize>>,
    >,
    rm_rbm: bool,
) -> Result<Vec<f64>, Box<dyn Error>> {
    let mut fem =
        FEM::from_pickle("data/20210802_0755_MT_mount_v202104_FSM/static_reduction_model.73.pkl")?;
    let n_io = (fem.n_inputs(), fem.n_outputs());
    //    println!("{}", fem);
    fem.keep_inputs(&[sid - 1]);
    fem.keep_outputs_by(&[sid, 25], |x| {
        x.descriptions.contains(&format!("M1-S{}", sid))
    });

    let nodes = fem.outputs[sid]
        .as_ref()
        .unwrap()
        .get_by(|x| x.properties.location.as_ref().map(|x| x.to_vec()))
        .into_iter()
        .flatten()
        .collect::<Vec<f64>>();
    let n_nodes = nodes.len() / 3;

    let gain = fem.reduced_static_gain(n_io).unwrap();
    let buf = gain.as_slice().to_vec();
    let g = na::DMatrix::from_column_slice(gain.nrows(), gain.ncols(), &buf);
    let y = g * force;

    let (shape, rbm) = y.as_slice().split_at(n_nodes);
    if rm_rbm {
        let (t_xyz, r_xyz) = rbm.split_at(3);
        let rxyz_surface = {
            let q = Quaternion::unit(r_xyz[2], Vector::k())
                * Quaternion::unit(r_xyz[1], Vector::j())
                * Quaternion::unit(r_xyz[0], Vector::i());
            let trans_nodes: Vec<f64> = nodes
                .chunks(3)
                .flat_map(|x| {
                    let p: Quaternion = Vector::from(x).into();
                    let w: Quaternion = &q * p * &q.complex_conjugate();
                    Vec::<f64>::from(Vector::from(w.vector_as_slice()))
                })
                .collect();
            trans_nodes
                .chunks(3)
                .map(|x| x[2])
                .zip(nodes.chunks(3).map(|x| x[2]))
                .map(|(z_rbm, z)| z_rbm - z)
                .collect::<Vec<f64>>()
        };
        Ok(shape
            .iter()
            .zip(rxyz_surface.iter())
            .map(|(a, r)| a - r - t_xyz[2])
            .collect::<Vec<f64>>())
    } else {
        Ok(shape.to_owned())
    }
}

struct DataPoint {
    point: [f64; 2],
    data: f64,
}
impl HasPosition for DataPoint {
    type Point = [f64; 2];
    fn position(&self) -> [f64; 2] {
        self.point
    }
}

fn eigen_modes(
    sid: usize,
    rbm_eigen_forces: Option<&Matrix>,
    plot: Option<bool>,
) -> Result<U_BM2F, Box<dyn Error>> {
    let mut fem =
        FEM::from_pickle("data/20210802_0755_MT_mount_v202104_FSM/static_reduction_model.73.pkl")?;
    let n_io = (fem.n_inputs(), fem.n_outputs());
    //    println!("{}", fem);
    fem.keep_inputs(&[sid - 1]);
    fem.keep_outputs_by(&[sid, 25], |x| {
        x.descriptions.contains(&format!("M1-S{}", sid))
    });
    println!("{}", fem);

    let nodes = fem.outputs[sid]
        .as_ref()
        .unwrap()
        .get_by(|x| x.properties.location.as_ref().map(|x| x.to_vec()))
        .into_iter()
        .flatten()
        .collect::<Vec<f64>>();
    let n_nodes = nodes.len() / 3;

    let gain = fem.reduced_static_gain(n_io).unwrap();
    let mut m1s_influences = vec![];

    for col in gain.column_iter() {
        let (shape, rbm) = col.as_slice().split_at(n_nodes);
        let (t_xyz, r_xyz) = rbm.split_at(3);
        let rxyz_surface = {
            let q = Quaternion::unit(r_xyz[2], Vector::k())
                * Quaternion::unit(r_xyz[1], Vector::j())
                * Quaternion::unit(r_xyz[0], Vector::i());
            let trans_nodes: Vec<f64> = nodes
                .chunks(3)
                .flat_map(|x| {
                    let p: Quaternion = Vector::from(x).into();
                    let w: Quaternion = &q * p * &q.complex_conjugate();
                    Vec::<f64>::from(Vector::from(w.vector_as_slice()))
                })
                .collect();
            trans_nodes
                .chunks(3)
                .map(|x| x[2])
                .zip(nodes.chunks(3).map(|x| x[2]))
                .map(|(z_rbm, z)| z_rbm - z)
                .collect::<Vec<f64>>()
        };
        m1s_influences.extend(
            shape
                .iter()
                .zip(rxyz_surface.iter())
                .map(|(a, r)| a - r - t_xyz[2])
                .collect::<Vec<f64>>(),
        );
    }

    /*
       gain.as_slice()
           .chunks(n_nodes + 6)
           .enumerate()
           .for_each(|(k, col)| {
               heatmap(
                   &format!("data/m1_modes/zonal/m1_if-{:03}.svg", k + 1),
                   &col[..n_nodes],
                   &nodes,
               );
           });

       m1s_influences
           .as_slice()
           .chunks(n_nodes)
           .enumerate()
           .for_each(|(k, col)| {
               heatmap(
                   &format!("data/m1_modes/zonal_wo-rbm/m1_if-{:03}.svg", k + 1),
                   col,
                   &nodes,
               );
           });
    */
    let mat =
        na::DMatrix::from_column_slice(n_nodes, m1s_influences.len() / n_nodes, &m1s_influences);
    let m1s_svd = match rbm_eigen_forces {
        Some(v_rbm_t) => {
            let m1s_svd = mat.svd(true, true);
            let v = m1s_svd.v_t.as_ref().unwrap().transpose();
            let v_wo_rbm = &v - (v_rbm_t.transpose() * (v_rbm_t * &v));
            let mat = m1s_svd.u.as_ref().unwrap()
                * na::DMatrix::from_diagonal(&m1s_svd.singular_values)
                * v_wo_rbm.transpose();
            mat.svd(true, true)
        }
        None => mat.svd(true, true),
    };

    println!("U: {:?}", m1s_svd.u.as_ref().unwrap().shape());

    let inv_s = m1s_svd.singular_values.map(|x| x.recip());
    println!("S^-1: [{:?}]", inv_s.shape());
    println!(
        "V: [{:?}]",
        m1s_svd.v_t.as_ref().unwrap().transpose().shape()
    );
    let bm2force = m1s_svd.v_t.as_ref().unwrap().transpose() * na::DMatrix::from_diagonal(&inv_s);

    if plot.is_some() {
        m1s_svd
            .u
            .as_ref()
            .unwrap()
            .column_iter()
            .enumerate()
            .for_each(|(k, col)| {
                heatmap(
                    &format!("data/m1_modes/modal/m1_s{}_if-{:03}.svg", sid, k + 1),
                    col.as_slice(),
                    &nodes,
                    None,
                );
            });
    }

    Ok((m1s_svd.u.unwrap(), bm2force, nodes))
}

fn gridding(bm: &[f64], nodes: &[f64]) -> Result<Vec<f64>, Box<dyn Error>> {
    let n = 201;
    let width = 8.4;
    let delta = 8.4 / (n - 1) as f64;
    let n_node = nodes.len() / 3;
    let n_mode = bm.len() / n_node;
    let mut data_gridded: Vec<f64> = Vec::with_capacity(n * n * n_mode);
    let mut modes = bm.chunks(n_node);
    for k in 0..n_mode {
        let mut delaunay = FloatDelaunayTriangulation::with_walk_locate();
        nodes
            .chunks(3)
            .zip(modes.next().unwrap().iter())
            .for_each(|(node, &data)| {
                delaunay.insert(DataPoint {
                    point: [node[0], node[1]],
                    data,
                });
            });
        for i in 0..n {
            let x = i as f64 * delta - width * 0.5;
            for j in 0..n {
                let y = j as f64 * delta - width * 0.5;
                if let Some(interpolated) = delaunay.nn_interpolation(&[x, y], |dp| dp.data) {
                    data_gridded.push(interpolated)
                } else {
                    return Err("Interpolation failed.".into());
                };
            }
        }
    }
    Ok(data_gridded)
}

fn bm2ceo() -> Result<(), Box<dyn Error>> {
    let n_max = 329 * 201 * 201;
    let u_grid: Vec<f64> = (1..=7)
        .flat_map(|sid| {
            let filename = format!("m1s{}bm.pkl", sid);
            let mut buf: Vec<f64> = pkl::from_reader(File::open(filename).unwrap()).unwrap();
            if sid < 7 {
                buf.truncate(n_max)
            } else {
                let n_buf = buf.len();
                buf.extend_from_slice(&vec![0f64; n_max - n_buf])
            }
            buf
        })
        .collect();
    let mut file = File::create("m1_eigen_modes.ceo")?;
    file.write_all(&201_i32.to_ne_bytes());
    file.write_all(&8.4_f64.to_ne_bytes());
    file.write_all(&7_i32.to_ne_bytes());
    file.write_all(&329_i32.to_ne_bytes());
    let s2b: Vec<_> = (0i32..7i32).flat_map(|i| i.to_ne_bytes()).collect();
    file.write_all(&s2b);
    let m: Vec<_> = u_grid.iter().flat_map(|x| x.to_ne_bytes()).collect();
    file.write_all(&m);
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    //    let rbm_eigen_forces = m1_actuators_2_rbm(sid)?;
    bm2ceo()?;
    /*
        let (u_bm2f, u_grid): (Vec<_>, Vec<_>) = (1..=7)
            .into_par_iter()
            .map(|sid| {
                let u_bm2f = eigen_modes(sid, Some(&m1_actuators_2_rbm(sid).unwrap()), None).unwrap();
                let u_grid = gridding(&u_bm2f.0.as_slice(), &u_bm2f.2)
                    .map_err(|e| println!("{:?}", e))
                    .unwrap();
                (u_bm2f, u_grid)
            })
            .unzip();

        u_grid.iter().enumerate().for_each(|(k, u_grid)| {
            let filename = format!("m1s{}bm.pkl", k + 1);
            let mut file = File::create(filename).unwrap();
            pkl::to_writer(&mut file, u_grid, true).unwrap();
        });

        let sid = 2;
        let u_bm2f_sid = &u_bm2f[sid - 1];
        let b: na::DVector<f64> = na::DVector::from_iterator(
            u_bm2f_sid.0.ncols(),
            (0..u_bm2f_sid.0.ncols()).map(|k| {
                if k < 42 {
                    let r: f64 = rand::random();
                    (2f64 * r - 1f64) * 1e-6
                } else {
                    0f64
                }
            }),
        );
        let surface = &u_bm2f_sid.0 * &b;
        heatmap(
            "surf_from_u.svg",
            surface.as_slice(),
            &u_bm2f_sid.2,
            Some(1e6),
        );
        println!("bm2f shape: {:?}", u_bm2f_sid.1.shape());
        let force = &u_bm2f_sid.1 * b;
        let fem_surface = apply_force(sid, force, false)?;
        heatmap(
            "surf_from_fem.svg",
            fem_surface.as_slice(),
            &u_bm2f_sid.2,
            Some(1e6),
        );

        let surface_error_rms = (surface
            .iter()
            .zip(fem_surface.iter())
            .map(|(a, b)| a - b)
            .fold(0f64, |s, x| s + x * x)
            / surface.len() as f64)
            .sqrt();
        println!("WFE RMS: {:8.6}nm", surface_error_rms * 1e9);

        let (mode, coefs, _) = zernike::filter(
            surface.as_slice(),
            u_bm2f_sid.2.chunks(3).map(|xyz| (xyz[0], xyz[1])),
            4,
        );
        println!(
            "Zern. coefs. : {:#?}",
            coefs.iter().map(|x| x * 1e6).collect::<Vec<f64>>()
        );
        let (mode, coefs, _) = zernike::filter(
            fem_surface.as_slice(),
            u_bm2f_sid.2.chunks(3).map(|xyz| (xyz[0], xyz[1])),
            4,
        );
        println!(
            "Zern. coefs. : {:#?}",
            coefs.iter().map(|x| x * 1e6).collect::<Vec<f64>>()
        );
    */
    Ok(())
}
