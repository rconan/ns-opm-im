use super::Config;
use colorous;
use plotters::prelude::*;

pub struct Combo {}
pub enum Kind {
    Plot,
    Scatter,
}
pub type Complot = (
    Vec<Box<(dyn Iterator<Item = (f64, Vec<f64>)> + 'static)>>,
    Vec<Kind>,
    Option<Config<'static>>,
);
impl<'a> From<Complot> for Combo {
    fn from((iters, draws, config): Complot) -> Self {
        let config = config.unwrap_or_default();
        let filename = config
            .filename
            .unwrap_or_else(|| "complot-plot.svg".to_string());

        let fig = SVGBackend::new(&filename, (768, 512)).into_drawing_area();
        fig.fill(&WHITE).unwrap();

        let xrange = config.xaxis.range.unwrap();
        let yrange = config.yaxis.range.unwrap();

        let mut chart = ChartBuilder::on(&fig)
            .set_label_area_size(LabelAreaPosition::Left, 50)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .margin(10)
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

        for (iter, draw) in iters.into_iter().zip(draws.into_iter()) {
            let mut colors = colorous::TABLEAU10.iter().cycle();
            let xy: Vec<_> = iter.collect();
            let n_y = xy.iter().nth(0).unwrap().1.len();
            let data: Vec<_> = xy
                .into_iter()
                .flat_map(|(x, y)| y.into_iter().map(|y| (x, y)).collect::<Vec<(f64, f64)>>())
                .collect();
            match draw {
                Kind::Scatter => {
                    for k in 0..n_y {
                        let this_color = colors.next().unwrap().as_tuple();
                        chart
                            .draw_series(data.iter().skip(k).step_by(n_y).cloned().map(|point| {
                                Circle::new(
                                    point,
                                    3,
                                    RGBColor(this_color.0, this_color.1, this_color.2).filled(),
                                )
                            }))
                            .unwrap();
                    }
                }
                Kind::Plot => {
                    for k in 0..n_y {
                        let this_color = colors.next().unwrap().as_tuple();
                        chart
                            .draw_series(LineSeries::new(
                                data.iter().skip(k).step_by(n_y).cloned(),
                                RGBColor(this_color.0, this_color.1, this_color.2),
                            ))
                            .unwrap();
                    }
                }
            }
        }
        Combo {}
    }
}
