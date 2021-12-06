use complot::{Axis, Config, Plot};
use crseo::{from_opticals, OpticalSensitivities};
use dosio::{ios, Dos, IOVec, IO};
use fem::{
    dos::{DiscreteModalSolver, DiscreteStateSpace, Exponential},
    FEM,
};
use geotrans::{Quaternion, Vector};
use indicatif::{ProgressBar, ProgressStyle};
use log;
use m1_ctrl as m1;
use mount_ctrl as mount;
use nalgebra as na;
use ns_opm_im::control::{Control, Delay};
use serde_pickle as pickle;
use skyangle::Conversion;
use std::{fs::File, io::BufWriter, path::Path, time::Instant};
use windloading::WindLoads;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    let sim_duration = 30f64;
    let sampling_rate = 1e3;
    let wfs_delay = 10f64;
    let wfs_sample_delay = (sampling_rate * wfs_delay) as usize;
    let aco_delay = 1f64;
    let aco_sample_delay = (sampling_rate * aco_delay) as usize;
    let tiptilt_exposure_time = 5e-3;
    let tiptilt_sample_rate = (tiptilt_exposure_time * sampling_rate) as usize;
    let aco_exposure_time = 30f64;
    let aco_sample_rate = (aco_exposure_time * sampling_rate) as usize;
    let m1_sampling_rate = 100f64;
    let m1_sample_rate = (sampling_rate / m1_sampling_rate) as usize;
    println!(
        r##"
NS OPM IM Timing:
 - duration                        : {:}s,
 - sampling rate                   : {:}Hz,
 - Tip-tilt loop delay             : {:}s ({:} sample),
 - Tip-tilt loop sampling rate     : {:}Hz ({:} sample),
 - Active Optics loop delay        : {:}s ({:} sample),
 - Active Optics loop sampling rate: {:}Hz ({:} sample),
 - M1 force loop sampling rate     : {:}Hz ({:} sample),
"##,
        sim_duration,
        sampling_rate,
        wfs_delay,
        wfs_sample_delay,
        tiptilt_exposure_time.recip(),
        tiptilt_sample_rate,
        aco_delay + wfs_delay,
        aco_sample_delay + wfs_sample_delay,
        aco_exposure_time.recip(),
        aco_sample_rate,
        m1_sampling_rate,
        m1_sample_rate
    );

    // WIND LOADS
    println!("Loading wind loads ...");
    let now = Instant::now();
    let mut wind_loading = WindLoads::from_pickle(
        //        Path::new("data").join("b2019_0z_30az_os_7ms.wind_loads_1kHz_100-400s.pkl"),
        Path::new("data")
            .join("windloading")
            .join("zen30az000_OS7")
            .join("zen30az000_OS7.wind_loads_1kHz_000-400s.pkl"),
    )?
    .range(0.0, sim_duration)
    //    .cring()?
    //    .gir()?
    //    .truss()?
    .m1_segments()?
    //.m1_cell()?
    //    .m2_asm_topend()?
    //    .m2_segments_into(ios!(MCM2LclForce6F))?
    .build()?;
    println!("... in {}ms", now.elapsed().as_millis());

    // MOUNT
    let mut mnt_drives = mount::drives::Controller::new();
    let mut mnt_ctrl = mount::controller::Controller::new();
    // M1
    let mut m1_hardpoints = m1::hp_dynamics::Controller::new();
    let mut m1_load_cells = m1::hp_load_cells::Controller::new();
    let mut m1_actuators = m1::actuators::M1ForceLoops::new();
    // M2
    let mut fsm_positionner = fsm::positionner::Controller::new();
    let mut fsm_piezostack = fsm::piezostack::Controller::new();
    // Optical System
    // - 1 sample delay
    // let mut delay_1 = Delay::new(1, 14);
    // - tipt-tilt
    // let mut os_tiptilt = fsm::tiptilt::Controller::new();
    //FEM
    let now = Instant::now();
    let (m1_segments_surface_nodes, mut fem): (Vec<Vec<f64>>, DiscreteModalSolver<Exponential>) = {
        let fem = FEM::from_env()?;
        //println!("FEM:\n{}", fem);
        (
            fem.outputs
                .iter()
                .skip(1)
                .step_by(3)
                .take(7)
                .filter_map(|io| io.as_ref())
                .map(|io| {
                    io.get_by(|x| x.properties.location.as_ref().map(|x| x.to_vec()))
                        .into_iter()
                        .flatten()
                        .collect::<Vec<f64>>()
                })
                .collect(),
            DiscreteStateSpace::from(fem)
                .sampling(sampling_rate)
                .proportional_damping(2. / 100.)
                .max_eigen_frequency(75f64)
                //.truncate_hankel_singular_values(1e-5)
                .inputs_from(&[
                    &wind_loading,
                    &mnt_drives,
                    &m1_hardpoints,
                    &m1_actuators,
                    &fsm_positionner,
                    &fsm_piezostack,
                ])
                .outputs_to(&[&mnt_ctrl])
                .outputs(vec![ios!(OSSHardpointD)])
                .outputs(ios!(MCM2SmHexD, MCM2PZTD))
                .outputs(ios!(OSSM1Lcl, MCM2Lcl6D))
                .build()?,
        )
    };
    println!("... in {}ms", now.elapsed().as_millis());
    println!("{}", fem);

    let optical_sensitivities = OpticalSensitivities::load()?;
    let rxy_2_stt = &optical_sensitivities[3].m2_rxy()?;
    let opticals = from_opticals(&optical_sensitivities[1..4]);

    // I/O initialization
    let mut m1_actuators_forces = Some(ios!(
        M1ActuatorsSegment1(vec![0f64; 335]),
        M1ActuatorsSegment2(vec![0f64; 335]),
        M1ActuatorsSegment3(vec![0f64; 335]),
        M1ActuatorsSegment4(vec![0f64; 335]),
        M1ActuatorsSegment5(vec![0f64; 335]),
        M1ActuatorsSegment6(vec![0f64; 335]),
        M1ActuatorsSegment7(vec![0f64; 306])
    ));
    let m1_bending_modes = ios!(
        M1S1BMcmd(vec![0f64; 335]),
        M1S2BMcmd(vec![0f64; 335]),
        M1S3BMcmd(vec![0f64; 335]),
        M1S4BMcmd(vec![0f64; 335]),
        M1S5BMcmd(vec![0f64; 335]),
        M1S6BMcmd(vec![0f64; 335]),
        M1S7BMcmd(vec![0f64; 306])
    );
    let mut m1_hardpoints_forces = Some(vec![ios!(OSSHarpointDeltaF(vec![0f64; 42]))]);
    let m1_rbm_cmd: Option<Vec<IO<Vec<f64>>>> = None;
    let mut m2_pos_cmd = vec![ios!(M2poscmd(vec![0f64; 42]))];
    let mut pzt_cmd = Some(vec![ios!(TTcmd(vec![0f64; 21]))]);
    let mut fem_outputs = fem.zeroed_outputs();

    let n_step = wind_loading.n_sample;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    //    let mut m1_bm_logs = Vec::<Vec<f64>>::with_capacity(n_step * (329 * 6 + 300));
    let mut m2_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    //    let mut wfs_log = Vec::<f64>::with_capacity(n_step * 14 / wfs_sample_rate);
    let mut data_log = Vec::<f64>::with_capacity(n_step * 23);
    let pb = ProgressBar::new(n_step as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{duration_precise}] {bar:60.cyan/blue} [{eta_precise}]")
            .progress_chars("=|-"),
    );
    let mut k = 0;
    let now = Instant::now();
    while let Some(mut fem_forces) = wind_loading.outputs() {
        //for k in 0..n_step {
        pb.inc(1);
        //let mut fem_forces = vec![];
        // MOUNT
        fem_outputs
            .pop_these(ios!(
                OSSElEncoderAngle,
                OSSAzEncoderAngle,
                OSSRotEncoderAngle
            ))
            .and_then(|mut mnt_encdr| {
                mnt_ctrl
                    .in_step_out(Some(mnt_encdr.clone()))
                    .unwrap()
                    .and_then(|mut mnt_cmd| {
                        mnt_cmd.append(&mut mnt_encdr);
                        mnt_drives.in_step_out(Some(mnt_cmd)).unwrap()
                    })
            })
            .map(|mut mount_drives_forces| {
                fem_forces.append(&mut mount_drives_forces);
            });
        // M1
        if k % m1_sample_rate == 0 {
            m1_actuators_forces = fem_outputs
                .pop_these(vec![ios!(OSSHardpointD)])
                .and_then(|mut hp_d| {
                    m1_hardpoints_forces.as_mut().and_then(|hp_f| {
                        hp_d.extend_from_slice(hp_f);
                        m1_load_cells.in_step_out(Some(hp_d)).unwrap()
                    })
                })
                .and_then(|mut hp_lc| {
                    hp_lc.extend_from_slice(&m1_bending_modes);
                    m1_actuators.in_step_out(Some(hp_lc)).unwrap()
                })
        };
        m1_actuators_forces.as_mut().map(|x| {
            fem_forces.extend_from_slice(x);
        });
        if let Some(ref value) = m1_rbm_cmd {
            m1_hardpoints_forces = m1_hardpoints.in_step_out(Some(value.to_owned())).unwrap();
        }
        m1_hardpoints_forces.as_mut().map(|x| {
            fem_forces.extend_from_slice(x);
        });
        // M2
        //  - positioner
        fem_outputs
            .pop_these(vec![ios!(MCM2SmHexD)])
            .and_then(|mut hex_d| {
                hex_d.extend_from_slice(&mut m2_pos_cmd);
                fsm_positionner.in_step_out(Some(hex_d)).unwrap()
            })
            .map(|mut fsm_positionner_forces| {
                fem_forces.append(&mut fsm_positionner_forces);
            });
        //  - piezostack
        fem_outputs
            .pop_these(vec![ios!(MCM2PZTD)])
            .and_then(|mut pzt_d| {
                pzt_cmd
                    .as_mut()
                    .map(|pzt_cmd| pzt_d.extend_from_slice(pzt_cmd));
                fsm_piezostack.in_step_out(Some(pzt_d)).unwrap()
            })
            .map(|mut fsm_piezostack_forces| {
                fem_forces.append(&mut fsm_piezostack_forces);
            });
        // FEM
        fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap().unwrap();
        let m1_rbm: Option<Vec<f64>> = fem_outputs[ios!(OSSM1Lcl)].clone().into();

        k += 1;
        // On-axis optical sensitivities
        let tspst = &opticals
            * na::DVector::from_iterator(
                84,
                Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)])
                    .into_iter()
                    .chain(Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).into_iter())
                    .flatten()
                    .into_iter(),
            );
        data_log.extend(&tspst);
        // LOGS
        Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)]).map(|x| m1_logs.push(x));
        Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).map(|x| m2_logs.push(x));
    }
    pb.finish();
    let eta = now.elapsed();
    println!(
        "{} sample wind loads played in [{}] ({:.3}ms/step)",
        n_step,
        humantime::format_duration(eta),
        eta.as_millis() as f64 / n_step as f64
    );

    println!("file: {:?}", Path::new(file!()).file_name());
    let root = "windloading";

    let (tilts, segments): (Vec<_>, Vec<_>) = data_log
        .iter()
        .cloned()
        .collect::<Vec<f64>>()
        .chunks(23)
        .map(|data| {
            let (a, b) = data.split_at(2);
            (a.to_vec(), b.to_vec())
        })
        .unzip();
    let (segments_piston, _segments_tilts): (Vec<_>, Vec<_>) = segments
        .iter()
        .flatten()
        .cloned()
        .collect::<Vec<f64>>()
        .chunks(21)
        .map(|data| {
            let (a, b) = data.split_at(7);
            (a.to_vec(), b.to_vec())
        })
        .unzip();

    let n_skip = 5_000;
    Plot::from((
        tilts
            .into_iter()
            .skip(n_skip)
            .flatten()
            .collect::<Vec<f64>>()
            .chunks(2)
            .enumerate()
            .map(|(k, tilts)| {
                (
                    (k + n_skip) as f64 / sampling_rate,
                    tilts.iter().map(|x| x.to_mas()).collect::<Vec<f64>>(),
                )
            }),
        Some(
            Config::new()
                .filename(format!("{}_tilts.svg", root))
                .xaxis(Axis::new().label("Time [s]"))
                .yaxis(Axis::new().label("Tip-tilt [mas]")),
        ),
    ));

    Plot::from((
        segments_piston
            .into_iter()
            .skip(n_skip)
            .flatten()
            .collect::<Vec<f64>>()
            .chunks(7)
            .enumerate()
            .map(|(k, pistons)| {
                (
                    (k + n_skip) as f64 / sampling_rate,
                    pistons.iter().map(|x| x * 1e6).collect::<Vec<f64>>(),
                )
            }),
        Some(
            Config::new()
                .filename(format!("{}_pistons.svg", root))
                .xaxis(Axis::new().label("Time [s]"))
                .yaxis(Axis::new().label("Piston [micron]")),
        ),
    ));

    let mut m1_segment_rbm = vec![vec![]; 7];
    m1_logs.iter().for_each(|m1_rbms| {
        m1_rbms
            .chunks(6)
            .zip(&mut m1_segment_rbm)
            .for_each(|(c, r)| r.push(c))
    });
    m1_segment_rbm.into_iter().enumerate().for_each(|(sid, x)| {
        Plot::from((
            x.into_iter().enumerate().map(|(k, wfe_rms)| {
                (
                    k as f64 / sampling_rate,
                    wfe_rms.iter().map(|x| x * 1e6).collect::<Vec<f64>>(),
                )
            }),
            Some(
                Config::new()
                    .filename(&format!("{}_m1-s{}-rbm.svg", root, sid + 1))
                    .xaxis(Axis::new().label("Time [s]"))
                    .yaxis(Axis::new().label("WFE RMS [micron]")),
            ),
        ));
    });

    let mut m2_segment_rbm = vec![vec![]; 7];
    m2_logs.iter().for_each(|m2_rbms| {
        m2_rbms
            .chunks(6)
            .zip(&mut m2_segment_rbm)
            .for_each(|(c, r)| r.push(c))
    });
    m2_segment_rbm.into_iter().enumerate().for_each(|(sid, x)| {
        Plot::from((
            x.into_iter().enumerate().map(|(k, wfe_rms)| {
                (
                    k as f64 / sampling_rate,
                    wfe_rms.iter().map(|x| x * 1e6).collect::<Vec<f64>>(),
                )
            }),
            Some(
                Config::new()
                    .filename(&format!("{}_m2-s{}-rbm.svg", root, sid + 1))
                    .xaxis(Axis::new().label("Time [s]"))
                    .yaxis(Axis::new().label("WFE RMS [micron]")),
            ),
        ));
    });
    /*
        let file = File::create("ns-opm-im_gmt-logs.pkl")?;
        let mut buf = BufWriter::with_capacity(1_000_000, file);
        serde_pickle::to_writer(&mut buf, &(m1_logs, m2_logs, m1_bm_logs), true)?;
    */
    Ok(())
}
