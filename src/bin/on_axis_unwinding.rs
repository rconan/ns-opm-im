use complot::{Axis, Config, Plot};
use crseo::{dos::GmtOpticalModel, Builder, GMT, SOURCE};
use dosio::{ios, Dos, IOVec};
use indicatif::{ProgressBar, ProgressStyle};
use std::{fs::File, io::BufReader};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    let file = File::open("ns-opm-im_gmt-logs.pkl")?;
    let buf = BufReader::with_capacity(1_000_000, file);
    let (m1_rbm_logs, m2_rbm_logs, m1_bm_logs): (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>) =
        serde_pickle::from_reader(buf)?;

    let sim_duration = 60f64;
    let sampling_rate = 1e3;
    let m1_n_mode = 332;
    let soak_delta_temperature = 10f64;
    let n_px = 16 * 48 + 1;
    let mut gom = GmtOpticalModel::new()
        .gmt(GMT::new().m1(
            "m1_eigen-modes_raw-polishing_print-through_soak1deg",
            m1_n_mode,
        ))
        .source(SOURCE::new().pupil_sampling(n_px))
        .output(ios!(SrcSegmentWfeRms))
        .dome_seeing("b2019_0z_0az_os_7ms", sim_duration, sampling_rate, None)
        .await?
        .build()?;
    gom.gmt.a1 = (0..7)
        .flat_map(|_| {
            let mut a1 = vec![0f64; m1_n_mode];
            a1[m1_n_mode - 3] = 1f64;
            a1[m1_n_mode - 2] = 1f64;
            a1[m1_n_mode - 1] = soak_delta_temperature;
            a1
        })
        .collect();
    gom.gmt.reset();

    let n_skip = 0;
    let n_sample = m1_rbm_logs.len();
    let mut segment_wfe_rms = Vec::<Option<Vec<f64>>>::with_capacity(n_sample);
    let pb = ProgressBar::new((n_sample - n_skip) as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{duration_precise}] {bar:60.cyan/blue} [{eta_precise}]")
            .progress_chars("=|~"),
    );
    for ((m1_rbm, m2_rbm), m1_bm) in m1_rbm_logs
        .iter()
        .zip(&m2_rbm_logs)
        .zip(&m1_bm_logs)
        .skip(n_skip)
    {
        pb.inc(1);
        let mut data = gom
            .in_step_out(Some(ios!(
                OSSM1Lcl(m1_rbm.to_owned()),
                MCM2Lcl6D(m2_rbm.to_owned()),
                M1modes(m1_bm.to_owned())
            )))?
            .unwrap();
        segment_wfe_rms.push(data.pop_this(ios!(SrcSegmentWfeRms)).unwrap().into());
    }
    pb.finish();

    let sampling_rate = 1e3;
    Plot::from((
        segment_wfe_rms
            .into_iter()
            .filter_map(|x| x)
            .enumerate()
            .map(|(k, wfe_rms)| {
                (
                    (n_skip + k) as f64 / sampling_rate,
                    wfe_rms.iter().map(|x| x * 1e9).collect::<Vec<f64>>(),
                )
            }),
        Some(
            Config::new()
                .filename("ns-opm-im_segment-wfe-rms.svg")
                .xaxis(Axis::new().label("Time [s]"))
                .yaxis(Axis::new().label("Segment WFE RMS [nm]")),
        ),
    ));

    let phase: Vec<_> = gom
        .src
        .segment_phase()
        .iter()
        .map(|&x| x as f64 * 1e6)
        .collect();
    let _: complot::Heatmap = ((phase.as_slice(), (n_px, n_px)), None).into();

    Ok(())
}
