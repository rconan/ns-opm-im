use complot::{Axis, Config, Plot};
use crseo::{dos::GmtOpticalModel, Builder, GMT, SOURCE};
use dosio::{ios, Dos, IOVec};
use indicatif::{ProgressBar, ProgressStyle};
use skyangle::Conversion;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    let sim_duration = 310f64;
    let sampling_rate = 1e3;

    let atm_duration = 20f32;
    let atm_n_duration = Some((sim_duration / atm_duration as f64).ceil() as i32);
    let atm_sampling = 48 * 16 + 1;
    atm_n_duration.map(|atm_n_duration| {
        log::info!(
            "Atnmosphere duration: {}x{}={}s",
            atm_duration,
            atm_n_duration,
            atm_duration * atm_n_duration as f32,
        )
    });

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
        .atmosphere(
            crseo::ATMOSPHERE::new()
                .remove_turbulence_layer(0)
                .ray_tracing(
                    25.5,
                    atm_sampling,
                    20f32.from_arcmin() ,
                    atm_duration,
                    Some("atmosphere/ns-opm-im_atm.bin".to_string()),
                    atm_n_duration,
                ),
        )
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

    let os_gmt_state = Some(ios!(
        OSSM1Lcl(vec![0f64; 42]),
        MCM2Lcl6D(vec![0f64; 42]),
        M1modes(vec![0f64; m1_n_mode * 7])
    ));

    let n_skip = 0;
    let n_sample = (sampling_rate * sim_duration) as usize;
    let mut segment_wfe_rms = Vec::<Option<Vec<f64>>>::with_capacity(n_sample);
    let pb = ProgressBar::new((n_sample - n_skip) as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{duration_precise}] {bar:60.cyan/blue} [{eta_precise}]")
            .progress_chars("=|~"),
    );
    for k in 0..n_sample {
        pb.inc(1);
        if let Some(ref mut atm) = gom.atm {
            atm.secs = k as f64 / sampling_rate;
        }
        let mut data = gom.in_step_out(os_gmt_state.clone())?.unwrap();
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
                .filename("ns-opm-im_segment-wfe-rms_wo-aco.svg")
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
    let _: complot::Heatmap = (
        (phase.as_slice(), (n_px, n_px)),
        Some(complot::Config::new().filename("on_axis_wo-aco_wavefront.png")),
    )
        .into();

    Ok(())
}
