use fem::FEM;
use geotrans::{Quaternion, Vector};
use nalgebra as na;
use ns_opm_im::complot::{tri, Config};
use skyangle::Conversion;
use std::error::Error;
use triangle_rs as triangle;

fn heatmap(filename: &str, map: &[f64], nodes: &[f64]) {
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
        .zip(cells.into_iter().map(|x| x * 1e6));
    tri::Heatmap::from((data, Some(Config::new().filename(filename))));
}

fn main() -> Result<(), Box<dyn Error>> {
    let nodes = {
        let fem = FEM::from_pickle(
            "data/20210802_0755_MT_mount_v202104_FSM/static_reduction_model.73.pkl",
        )?;
        println!("{}", fem);
        fem.outputs[1]
            .as_ref()
            .unwrap()
            .get_by(|x| x.properties.location.as_ref().map(|x| x.to_vec()))
            .into_iter()
            .flatten()
            .collect::<Vec<f64>>()
    };
    let n_nodes = nodes.len() / 3;
    println!("nodes #: {}", n_nodes);

    let mut fem = FEM::from_pickle(
        "data/20210802_0755_MT_mount_v202104_FSM/xyz/static_reduction_model.73.pkl",
    )?;
    let n_io = (fem.n_inputs(), fem.n_outputs());
    println!("{}", fem);
    fem.keep_inputs(&[0]);
    let outputs_ids = [1usize, 2, 3, 39];
    fem.keep_outputs_by(&outputs_ids, |x| x.descriptions.contains("M1-S1"));
    println!("{}", fem);
    let n_outputs: Vec<usize> = outputs_ids
        .iter()
        .map(|&k| fem.outputs[k].as_ref().map(|o| o.len()).unwrap())
        .collect();
    println!("outputs size: {:?}", n_outputs);

    let fem_gain = fem.reduced_static_gain(n_io).unwrap();
    let u = na::DVector::from_element(fem.n_inputs(), 10f64);
    let y = &fem_gain * u;

    let y_range = |i: usize| {
        let start = n_outputs.iter().take(i).sum::<usize>();
        start..start + n_outputs[i]
    };
    let (t_xyz, r_xyz) = y.as_slice().get(y_range(3)).unwrap().split_at(3);
    println!(
        "Txyz [micron]: {:#?}",
        t_xyz.iter().map(|x| x * 1e6).collect::<Vec<f64>>()
    );
    println!(
        "Mean Tx [micron]: {:#?}",
        1e6 * y
            .as_slice()
            .get(y_range(1))
            .unwrap()
            .iter()
            .cloned()
            .sum::<f64>()
            / n_outputs[1] as f64
    );
    println!(
        "Mean Ty [micron]: {:#?}",
        1e6 * y
            .as_slice()
            .get(y_range(2))
            .unwrap()
            .iter()
            .cloned()
            .sum::<f64>()
            / n_outputs[2] as f64
    );
    println!(
        "Mean Tz [micron]: {:#?}",
        1e6 * y
            .as_slice()
            .get(y_range(0))
            .unwrap()
            .iter()
            .cloned()
            .sum::<f64>()
            / n_outputs[0] as f64
    );
    println!(
        "Rxyz [mas]: {:#?}",
        r_xyz.iter().map(|x| x.to_mas()).collect::<Vec<f64>>()
    );

    {
        let mode = y.as_slice().get(y_range(0)).unwrap();
        heatmap("m1-s1-surf.svg", &mode, &nodes);
    }
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
        let mode = trans_nodes
            .chunks(3)
            .map(|x| x[2])
            .zip(nodes.chunks(3).map(|x| x[2]))
            .map(|(z_rbm, z)| z_rbm - z)
            .collect::<Vec<f64>>();
        heatmap("m1-s1_rxyz-surf.svg", &mode, &nodes);
        mode
    };

    let no_rbm_surf = {
        let mode = y
            .as_slice()
            .get(y_range(0))
            .unwrap()
            .iter()
            .zip(rxyz_surface.iter())
            .map(|(a, r)| a - r - t_xyz[2])
            .collect::<Vec<f64>>();
        heatmap("m1-s1-surf-wo-rbm.svg", &mode, &nodes);
        mode
    };

    {
        let (mode, coefs, _) = zernike::filter(
            y.as_slice().get(y_range(0)).as_ref().unwrap(),
            nodes.chunks(3).map(|xyz| (xyz[0], xyz[1])),
            2,
        );
        println!(
            "Zern. coefs. : {:#?}",
            coefs.iter().map(|x| x * 1e6).collect::<Vec<f64>>()
        );
        heatmap("m1-s1-surf-wo-ptt.svg", &mode, &nodes);
        mode
    };

    {
        let (mode, coefs, _) =
            zernike::filter(&no_rbm_surf, nodes.chunks(3).map(|xyz| (xyz[0], xyz[1])), 4);
        println!(
            "Zern. coefs. : {:#?}",
            coefs.iter().map(|x| x * 1e6).collect::<Vec<f64>>()
        );
        heatmap("m1-s1-zern-surf.svg", &mode, &nodes);
    }

    Ok(())
}
