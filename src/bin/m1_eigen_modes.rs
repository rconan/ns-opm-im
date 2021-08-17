use fem::FEM;
use indicatif::ProgressBar;
use ns_opm_im::complot::{tri, Config};
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
    let mut fem =
        FEM::from_pickle("data/20210802_0755_MT_mount_v202104_FSM/static_reduction_model.73.pkl")?;
    let n_io = (fem.n_inputs(), fem.n_outputs());
    println!("{}", fem);
    fem.keep_inputs(&[0]);
    fem.keep_outputs_by(&[1], |x| x.descriptions.contains("M1-S1"));
    println!("{}", fem);

    let nodes = fem.outputs[1]
        .as_ref()
        .unwrap()
        .get_by(|x| x.properties.location.as_ref().map(|x| x.to_vec()))
        .into_iter()
        .flatten()
        .collect::<Vec<f64>>();
    let n_nodes = nodes.len() / 3;

    let m1s_influences = fem.reduced_static_gain(n_io).unwrap();
    heatmap(
        "m1_if.svg",
        m1s_influences.as_slice().chunks(n_nodes).nth(200).unwrap(),
        &nodes,
    );

    /*
        let m1s_svd = m1s_influences.svd(true, true);

        println!("U: {:?}", m1s_svd.u.as_ref().unwrap().shape());

        let nodes = fem.outputs[1]
            .as_ref()
            .unwrap()
            .get_by(|x| x.properties.location.as_ref().map(|x| x[..2].to_vec()))
            .into_iter()
            .flatten()
            .collect::<Vec<f64>>();

        let delaunay = triangle::Builder::new().set_tri_points(nodes).build();
        /*    tri::Mesh::from((
            delaunay.triangle_vertex_iter(),
            Some(Config::new().filename("m1-s1.svg")),
        ));*/

        let pb = ProgressBar::new(fem.n_inputs() as u64);
        for (k, mode) in m1s_svd
            .u
            .as_ref()
            .unwrap()
            .as_slice()
            .chunks(fem.n_outputs())
            .enumerate()
        {
            pb.inc(1);
            let cells: Vec<f64> = delaunay
                .triangle_iter()
                .map(|t| t.iter().fold(0., |a, &i| a + mode[i] / 3.))
                .collect();

            let data = delaunay.triangle_vertex_iter().zip(cells.into_iter());
            let filename = format!("data/m1-eigen-modes/s1/mode_{:03}.svg", k + 1);
            tri::Heatmap::from((data, Some(Config::new().filename(filename))));
        }
        pb.finish();
    */
    Ok(())
}
