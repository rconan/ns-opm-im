use super::Config;
use colorous;
use plotters::{coord::types::RangedCoordf64, prelude::*};

pub struct Scatter {}

impl<'a, I: Iterator<Item = (f64, Vec<f64>)>> From<(I, Option<Config<'a>>)> for Scatter {
    fn from((iter, config): (I, Option<Config>)) -> Self {
        let config = config.unwrap_or_default();
        let filename = config
            .filename
            .unwrap_or_else(|| "complot-plot.svg".to_string());

        let fig = SVGBackend::new(&filename, (768, 768)).into_drawing_area();
        fig.fill(&WHITE).unwrap();
        let xy: Vec<_> = iter.collect();

        let (x_max, y_max) = xy.iter().cloned().fold(
            (f64::NEG_INFINITY, f64::NEG_INFINITY),
            |(fx, fy), (x, y)| {
                (
                    fx.max(x),
                    fy.max(y.iter().cloned().fold(f64::NEG_INFINITY, |fy, y| fy.max(y))),
                )
            },
        );
        let (x_min, y_min) =
            xy.iter()
                .cloned()
                .fold((f64::INFINITY, f64::INFINITY), |(fx, fy), (x, y)| {
                    (
                        fx.min(x),
                        fy.min(y.iter().cloned().fold(f64::INFINITY, |fy, y| fy.min(y))),
                    )
                });

        let xrange = if let Some(xrange) = config.xaxis.range {
            xrange
        } else {
            x_min..x_max
        };
        let yrange = if let Some(yrange) = config.yaxis.range {
            yrange
        } else {
            y_min..y_max
        };

        let mut chart = ChartBuilder::on(&fig)
            //            .set_label_area_size(LabelAreaPosition::Left, 50)
            //            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .margin(20)
            .build_cartesian_2d(xrange, yrange)
            .unwrap();
        let mut mesh = chart.configure_mesh();
        if let Some(value) = config.xaxis.label {
            mesh.x_desc(value);
        }
        if let Some(value) = config.yaxis.label {
            mesh.y_desc(value);
        }
        mesh.draw().unwrap();

        let n_y = xy.iter().nth(0).unwrap().1.len();
        let data: Vec<_> = xy
            .into_iter()
            .flat_map(|(x, y)| y.into_iter().map(|y| (x, y)).collect::<Vec<(f64, f64)>>())
            .collect();
        let mut colors = colorous::TABLEAU10.iter().cycle();
        for k in 0..n_y {
            let this_color = colors.next().unwrap().as_tuple();
            chart
                .draw_series(data.iter().skip(k).step_by(n_y).cloned().map(|point| {
                    Circle::new(point, 3, RGBColor(this_color.0, this_color.1, this_color.2))
                }))
                .unwrap();
        }
        Scatter {}
    }
}
