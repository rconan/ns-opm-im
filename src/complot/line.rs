use super::{Config, Utils};
use colorous;
use plotters::prelude::*;
use std::iter::FromIterator;

pub struct Plot {}
impl Utils for Plot {}
/*impl<'a> FromIterator<f64> for Plot {
    fn from_iter<I: IntoIterator<Item = f64>>(iter: I) -> Self {
        let fig = SVGBackend::new("plot.svg", (768, 512)).into_drawing_area();
        fig.fill(&WHITE).unwrap();
        let xy: Vec<_> = iter
            .into_iter()
            .collect::<Vec<f64>>()
            .chunks(2)
            .map(|xy| (xy[0], xy[1]))
            .collect();
        let (x_max, y_max) = xy.iter().cloned().fold(
            (f64::NEG_INFINITY, f64::NEG_INFINITY),
            |(fx, fy), (x, y)| (fx.max(x), fy.max(y)),
        );
        let (x_min, y_min) = xy
            .iter()
            .cloned()
            .fold((f64::INFINITY, f64::INFINITY), |(fx, fy), (x, y)| {
                (fx.min(x), fy.min(y))
            });
        let mut ax = chart([x_min, x_max, y_min, y_max], &fig);
        ax.draw_series(xy.into_iter().map(|xy| Circle::new(xy, 3, RED.filled())))
            .unwrap();
        Plot {}
    }
}*/
impl<'a> FromIterator<(f64, Vec<f64>)> for Plot {
    fn from_iter<I: IntoIterator<Item = (f64, Vec<f64>)>>(iter: I) -> Self {
        let fig = SVGBackend::new("plot.svg", (768, 512)).into_drawing_area();
        fig.fill(&WHITE).unwrap();
        let xy: Vec<_> = iter.into_iter().collect();
        let (x_max, y_max) = Self::xy_max(&xy);
        let (x_min, y_min) = Self::xy_min(&xy);
        let mut chart = ChartBuilder::on(&fig)
            .set_label_area_size(LabelAreaPosition::Left, 50)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .margin(10)
            .build_cartesian_2d(x_min..x_max, y_min..y_max)
            .unwrap();
        let n_y = xy[0].1.len();
        let data: Vec<_> = xy
            .into_iter()
            .flat_map(|(x, y)| y.into_iter().map(|y| (x, y)).collect::<Vec<(f64, f64)>>())
            .collect();
        let mut colors = colorous::TABLEAU10.iter().cycle();
        for k in 0..n_y {
            let this_color = colors.next().unwrap().as_tuple();
            chart
                .draw_series(LineSeries::new(
                    data.iter().skip(k).step_by(n_y).cloned(),
                    RGBColor(this_color.0, this_color.1, this_color.2),
                ))
                .unwrap();
        }
        Plot {}
    }
}

impl<'a, I: Iterator<Item = (f64, Vec<f64>)>> From<(I, Option<Config<'a>>)> for Plot {
    fn from((iter, config): (I, Option<Config>)) -> Self {
        let mut config = config.unwrap_or_default();
        let filename = config
            .filename
            .unwrap_or_else(|| "complot-plot.svg".to_string());

        let fig = SVGBackend::new(&filename, (768, 512)).into_drawing_area();
        fig.fill(&WHITE).unwrap();
        let xy: Vec<_> = iter.collect();
        let (x_max, y_max) = Self::xy_max(&xy);
        let (x_min, y_min) = Self::xy_min(&xy);

        let mut chart = ChartBuilder::on(&fig)
            .set_label_area_size(LabelAreaPosition::Left, 50)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .margin(10)
            .build_cartesian_2d(x_min..x_max, y_min..y_max)
            .unwrap();
        let mut mesh = chart.configure_mesh();
        if let Some(value) = config.xaxis.label {
            mesh.x_desc(value);
        }
        if let Some(value) = config.yaxis.label {
            mesh.y_desc(value);
        }
        mesh.draw().unwrap();

        let n_y = xy[0].1.len();
        let data: Vec<_> = xy
            .into_iter()
            .flat_map(|(x, y)| y.into_iter().map(|y| (x, y)).collect::<Vec<(f64, f64)>>())
            .collect();
        let mut colors = colorous::TABLEAU10.iter().cycle();
        for k in 0..n_y {
            let this_color = colors.next().unwrap().as_tuple();
            chart
                .draw_series(LineSeries::new(
                    data.iter().skip(k).step_by(n_y).cloned(),
                    RGBColor(this_color.0, this_color.1, this_color.2),
                ))
                .unwrap();
        }
        Plot {}
    }
}
