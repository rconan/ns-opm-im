use fem::FEM;
use geotrans::m1_any_to_oss;
use plotters::prelude::*;
use std::{error::Error, ops::Deref};
use structopt::StructOpt;
use triangle_rs as mesh;

#[derive(Debug, StructOpt)]
struct Opt {
    /// Segment ID #
    sid: usize,
}

struct Segment<'a> {
    label: &'a str,
    id: usize,
    fem_id: usize,
}
enum Segments<'a> {
    Center(Segment<'a>),
    Outer(Segment<'a>),
}
impl<'a> Segments<'a> {
    fn center() -> Self {
        Segments::Center(Segment {
            label: "M1-S7",
            id: 7,
            fem_id: 6,
        })
    }
    fn outer(label: &'a str, id: usize, fem_id: usize) -> Self {
        Segments::Outer(Segment { label, id, fem_id })
    }
    fn one() -> Self {
        Segments::outer("M1-S1", 1, 0)
    }
    fn seven() -> Self {
        Segments::center()
    }
}
impl<'a> Deref for Segments<'a> {
    type Target = Segment<'a>;

    fn deref(&self) -> &Self::Target {
        use Segments::*;
        match self {
            Center(value) => value,
            Outer(value) => value,
        }
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let opt = Opt::from_args();
    let segment = match opt.sid {
        1 => Segments::one(),
        7 => Segments::seven(),
        _ => unimplemented!(),
    };

    let mut fem = FEM::from_env()?;
    println!("FEM:\n{}", &fem);
    fem.keep_inputs(&[segment.fem_id]);
    fem.keep_outputs_by(&[24], |x| {
        x.descriptions.starts_with(segment.label) && x.descriptions.ends_with("on mirror")
    });
    println!("FEM:\n{}", &fem);

    let mut actuators_coordinates = fem.inputs[segment.fem_id]
        .as_ref()
        .unwrap()
        .get_by(|x| x.properties.location.as_ref().map(|x| x[..2].to_vec()));
    actuators_coordinates.dedup_by(|a, b| a[0] == b[0] && a[1] == b[1]);
    println!(
        "actuators_coordinates: {}/{}",
        actuators_coordinates.len(),
        actuators_coordinates[0].len()
    );

    let mut hardpoints_coordinates = fem.outputs[24]
        .as_ref()
        .unwrap()
        .get_by(|x| x.properties.location.as_ref().map(|x| x.to_vec()));
    //    println!("hardpoints_coordinates: {:#?}", hardpoints_coordinates);
    hardpoints_coordinates.iter_mut().for_each(|c| {
        let v = m1_any_to_oss(segment.id, c.as_slice());
        c[0] = (*v)[0];
        c[1] = (*v)[1];
        c[2] = (*v)[2];
    });
    //    println!("hardpoints_coordinates: {:#?}", hardpoints_coordinates);

    let mut nodes = actuators_coordinates
        .iter()
        .flatten()
        .cloned()
        .collect::<Vec<f64>>();
    nodes.extend(
        hardpoints_coordinates
            .iter()
            .flat_map(|xyz| xyz[..2].to_vec()),
    );
    let mut builder = mesh::Builder::new().set_tri_points(nodes);
    let rim_diameter = 8.365;
    let delta_rim = 30e-2;
    let n_rim = (std::f64::consts::PI * rim_diameter / delta_rim).round();
    let outer_rim: Vec<_> = (0..n_rim as usize)
        .flat_map(|i| {
            let o = 2. * std::f64::consts::PI * i as f64 / n_rim;
            let (s, c) = o.sin_cos();
            let radius = 0.5 * rim_diameter;
            vec![radius * c, radius * s]
        })
        .collect();
    builder.add_polygon(&outer_rim);
    if let Segments::Center(_) = segment {
        let rim_diameter = 2.39;
        let n_rim = (std::f64::consts::PI * rim_diameter / delta_rim).round();
        let inner_rim: Vec<_> = (0..n_rim as usize)
            .flat_map(|i| {
                let o = 2. * std::f64::consts::PI * i as f64 / n_rim;
                let (s, c) = o.sin_cos();
                let radius = 0.5 * rim_diameter;
                vec![radius * c, radius * s]
            })
            .collect();
        builder.add_polygon(&inner_rim).add_holes(0., 0.);
    }
    let delaunay = builder.set_switches("pDqa0.075").build();
    println!("{}", delaunay);

    /*    let fig = complot::canvas("m1-s1_delaunay.svg");
    let lim = 4.5_f64;
    let mut ax = complot::chart([-lim, lim, -lim, lim], &fig);
    delaunay.mesh(&delaunay.x(), &delaunay.y(), [0; 3], &mut ax);
     */
    let pkl_file = format!("{}_delaunay.pkl", segment.label.to_lowercase());
    delaunay.dump(&pkl_file)?;

    let fig_file = format!("{}_delaunay.svg", segment.label.to_lowercase());
    let fig = SVGBackend::new(&fig_file, (768, 768)).into_drawing_area();
    fig.fill(&WHITE).unwrap();

    let xyrange = -4.3..4.3;

    let mut chart = ChartBuilder::on(&fig)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .margin(20)
        .build_cartesian_2d(xyrange.clone(), xyrange)
        .unwrap();
    let mut mesh = chart.configure_mesh();
    mesh.draw().unwrap();

    delaunay
        .triangle_iter()
        .map(|t| {
            t.iter()
                .map(|&i| (delaunay.x()[i], delaunay.y()[i]))
                .collect::<Vec<(f64, f64)>>()
        })
        .into_iter()
        .for_each(|v| {
            chart
                .draw_series(LineSeries::new(
                    v.iter().cycle().take(4).map(|(x, y)| (*x, *y)),
                    &BLACK,
                ))
                .unwrap();
        });

    let mut colors = colorous::TABLEAU10.iter().cycle();
    let iter = actuators_coordinates
        .into_iter()
        .map(|xy| (xy[0], vec![xy[1]]));
    let this_color = colors.next().unwrap().as_tuple();
    let xy: Vec<_> = iter.collect();
    let n_y = xy.iter().nth(0).unwrap().1.len();
    let data: Vec<_> = xy
        .into_iter()
        .flat_map(|(x, y)| y.into_iter().map(|y| (x, y)).collect::<Vec<(f64, f64)>>())
        .collect();
    for k in 0..n_y {
        chart
            .draw_series(data.iter().skip(k).step_by(n_y).cloned().map(|point| {
                Circle::new(
                    point,
                    5,
                    RGBColor(this_color.0, this_color.1, this_color.2).filled(),
                )
            }))
            .unwrap();
    }
    let iter = hardpoints_coordinates
        .into_iter()
        .map(|xyz| (xyz[0], vec![xyz[1]]));
    let this_color = colors.next().unwrap().as_tuple();
    let xy: Vec<_> = iter.collect();
    let n_y = xy.iter().nth(0).unwrap().1.len();
    let data: Vec<_> = xy
        .into_iter()
        .flat_map(|(x, y)| y.into_iter().map(|y| (x, y)).collect::<Vec<(f64, f64)>>())
        .collect();
    for k in 0..n_y {
        chart
            .draw_series(data.iter().skip(k).step_by(n_y).cloned().map(|point| {
                Circle::new(
                    point,
                    10,
                    RGBColor(this_color.0, this_color.1, this_color.2).filled(),
                )
            }))
            .unwrap();
    }

    /*
        <Scatters as From<(
            Vec<Box<(dyn Iterator<Item = (f64, Vec<f64>)> + 'static)>>,
            Option<Config<'static>>,
        )>>::from((
            vec![
                Box::new(
                    actuators_coordinates
                        .into_iter()
                        .map(|xy| (xy[0], vec![xy[1]])),
                ),
                Box::new(
                    hardpoints_coordinates
                        .into_iter()
                        .map(|xyz| (xyz[0], vec![xyz[1]])),
                ),
            ],
            Some(
                Config::new()
                    .filename("m1-meshing.svg")
                    .axes(Axis::new().range(-4.3..4.3)),
            ),
        ));
    */
    Ok(())
}

/*
M1-S7
hardpoints_coordinates: [
    [
        1.57483,
        -2.66543,
        3.455,
    ],
    [
        -1.574825,
        -2.66543,
        3.455,
    ],
    [
        -3.095739,
        -0.031125999999999987,
        3.455,
    ],
    [
        -1.520913,
        2.6965509999999995,
        3.4550000000000005,
    ],
    [
        1.52091,
        2.6965509999999995,
        3.4550000000000005,
    ],
    [
        3.09574,
        -0.031125999999999987,
        3.455,
    ],
]
*/
