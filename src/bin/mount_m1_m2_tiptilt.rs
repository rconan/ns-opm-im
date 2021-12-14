use complot::{Axis, Config, Plot};
use crseo::{
    calibrations,
    dos::GmtOpticalSensorModel,
    from_opticals,
    shackhartmann::{Diffractive, Geometric},
    Builder, Calibration, OpticalSensitivities, ShackHartmann, SH48, SOURCE,
};
use dosio::{ios, Dos, IOVec};
use fem::{
    dos::{DiscreteStateSpace, Exponential},
    FEM,
};
use indicatif::{ProgressBar, ProgressStyle};
use m1_ctrl as m1;
use mount_ctrl as mount;
use nalgebra as na;
use simple_logger::SimpleLogger;
use skyangle::Conversion;
use std::{path::Path, time::Instant};
use windloading::WindLoads;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    SimpleLogger::new().env().init()?;

    // WIND LOADS
    println!("Loading wind loads ...");
    let mut wind_loading = WindLoads::from_pickle(
        Path::new("data").join("b2019_0z_30az_os_7ms.wind_loads_1kHz_100-400s.pkl"),
    )?
    .range(0.0, 15.0)
    .truss()?
    .build()?;

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
    // - tipt-tilt
    let mut os_tiptilt = fsm::tiptilt::Controller::new();
    //FEM
    let sampling_rate = 1e3;
    let mut fem = {
        let fem = FEM::from_env()?;
        println!("FEM:\n{}", fem);
        DiscreteStateSpace::from(fem)
    }
    .sampling(sampling_rate)
    .proportional_damping(2. / 100.)
    .inputs_from(&[
        &wind_loading,
        &mnt_drives,
        &m1_hardpoints,
        &fsm_positionner,
        &fsm_piezostack,
    ])
    .inputs(ios!(
        M1ActuatorsSegment1,
        M1ActuatorsSegment2,
        M1ActuatorsSegment3,
        M1ActuatorsSegment4,
        M1ActuatorsSegment5,
        M1ActuatorsSegment6,
        M1ActuatorsSegment7
    ))
    .outputs(ios!(OSSM1Lcl, MCM2Lcl6D))
    .outputs(ios!(
        OSSAzEncoderAngle,
        OSSElEncoderAngle,
        OSSRotEncoderAngle
    ))
    .outputs(vec![ios!(OSSHardpointD)])
    .outputs(ios!(MCM2SmHexD, MCM2PZTD))
    .build::<Exponential>()?;

    // CEO WFSS
    let exposure_time = 5e-3;
    let n_sensor = 1;
    let wfs_sample_rate = (exposure_time * sampling_rate) as usize;
    let wfs_delay = sampling_rate as usize * 10;
    type WfsType = Diffractive;
    let mut gosm = GmtOpticalSensorModel::<ShackHartmann<WfsType>, SH48<WfsType>>::new(Some(
        SOURCE::new().size(n_sensor).fwhm(6f64),
    ))
    .sensor(SH48::<WfsType>::new().n_sensor(n_sensor))
    /*        .atmosphere(crseo::ATMOSPHERE::new().ray_tracing(
        26.,
        520,
        0.,
        25.,
        Some("ns-opm-im_atm.bin".to_string()),
        Some(8),
    ))*/
    .build()?;

    println!("M1 mode: {}", gosm.gmt.get_m1_mode_type());
    println!("M2 mode: {}", gosm.gmt.get_m2_mode_type());
    println!("GS band: {}", gosm.src.get_photometric_band());

    let mut gmt2wfs = Calibration::new(
        &gosm.gmt,
        &gosm.src,
        SH48::<Geometric>::new().n_sensor(n_sensor),
    );
    let mirror = vec![calibrations::Mirror::M2];
    let segments = vec![vec![calibrations::Segment::Rxyz(1e-6, Some(0..2))]; 7];
    let now = Instant::now();
    gmt2wfs.calibrate(
        mirror,
        segments,
        calibrations::ValidLensletCriteria::OtherSensor(&mut gosm.sensor),
    );
    println!(
        "GTM 2 WFS calibration [{}x{}] in {}s",
        gmt2wfs.n_data,
        gmt2wfs.n_mode,
        now.elapsed().as_secs()
    );
    let poke_sum = gmt2wfs.poke.from_dev().iter().sum::<f32>();
    println!("Poke sum: {}", poke_sum);
    let wfs_2_m2rxy = gmt2wfs.qr();
    let mut m2_seg_rbm = vec![vec![0f64; 6]; 7];
    m2_seg_rbm[1][3] = 1e-6;
    m2_seg_rbm[4][4] = 1e-6;
    m2_seg_rbm[6][3] = 1e-6;
    m2_seg_rbm[6][4] = 1e-6;

    let m2_rbm = ios!(MCM2Lcl6D(m2_seg_rbm.into_iter().flatten().collect()));
    //    gosm.inputs(vec![m2_rbm.clone()]).unwrap().step();
    let y = gosm
        .in_step_out(Some(vec![m2_rbm.clone()]))?
        .and_then(|x| Into::<Option<Vec<f64>>>::into(x[0].clone()))
        .map(|x| x.into_iter().map(|x| x as f32).collect::<Vec<f32>>())
        .unwrap();

    let a = wfs_2_m2rxy.solve(&mut y.into());
    Vec::<f32>::from(a)
        .into_iter()
        .map(|x| x * 1e6)
        .collect::<Vec<f32>>()
        .chunks(2)
        .enumerate()
        .for_each(|x| println!("#{}: [{:+0.1},{:+0.1}]", 1 + x.0, x.1[0], x.1[1]));

    // I/O initialization
    let mut mount_drives_forces = Some(ios!(
        OSSAzDriveTorque(vec![0f64; 12]),
        OSSElDriveTorque(vec![0f64; 4]),
        OSSRotDriveTorque(vec![0f64; 4])
    ));
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
    let mut fsm_positionner_forces = Some(vec![ios!(MCM2SmHexF(vec![0f64; 84]))]);
    let mut fsm_piezostack_forces = Some(vec![ios!(MCM2PZTF(vec![0f64; 42]))]);

    let mut fem_forces = vec![ios!(OSSTruss6F(vec![0f64; 36]))];
    mount_drives_forces.as_mut().map(|x| {
        fem_forces.append(x);
    });
    m1_hardpoints_forces.as_mut().map(|x| {
        fem_forces.extend_from_slice(x);
    });
    m1_actuators_forces.as_mut().map(|x| {
        fem_forces.extend_from_slice(x);
    });
    fsm_positionner_forces.as_mut().map(|x| {
        fem_forces.append(x);
    });
    fsm_piezostack_forces.as_mut().map(|x| {
        fem_forces.append(x);
    });

    let optical_sensitivities = OpticalSensitivities::load()?;
    let rxy_2_stt = &optical_sensitivities[3].m2_rxy()?;
    let opticals = from_opticals(&optical_sensitivities[1..4]);

    let mut pzt_cmd = Some(vec![ios!(TTcmd(vec![0f64; 21]))]);

    let mut fem_outputs = fem.in_step_out(Some(fem_forces))?.unwrap();

    let n_step = wind_loading.n_sample;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m2_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    //    let mut wfs_log = Vec::<f64>::with_capacity(n_step * 14 / wfs_sample_rate);
    let mut data_log = Vec::<f64>::with_capacity(n_step * 23);
    let pb = ProgressBar::new(n_step as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{duration_precise}] {bar:60.cyan/blue} {pos:>7}/{len:7}")
            .progress_chars("=|-"),
    );
    let mut k = 0;
    let now = Instant::now();
    while let Some(mut fem_forces) = wind_loading.outputs() {
        //for k in 0..n_step {
        pb.inc(1);
        //    let mut fem_forces = Vec::<IO<Vec<f64>>>::with_capacity(n_fem_forces);
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
        if k % 10 == 0 {
            m1_actuators_forces = fem_outputs
                .pop_these(vec![ios!(OSSHardpointD)])
                .and_then(|mut hp_d| {
                    m1_hardpoints_forces.as_mut().and_then(|hp_f| {
                        hp_d.append(hp_f);
                        m1_load_cells.in_step_out(Some(hp_d)).unwrap()
                    })
                })
                .and_then(|mut hp_lc| {
                    hp_lc.append(&mut m1_bending_modes.clone());
                    m1_actuators.in_step_out(Some(hp_lc)).unwrap()
                })
        };
        m1_actuators_forces.as_mut().map(|x| {
            fem_forces.extend_from_slice(x);
        });
        m1_hardpoints_forces = m1_hardpoints
            .in_step_out(Some(vec![ios!(M1RBMcmd(vec![0f64; 42]))]))
            .unwrap();
        m1_hardpoints_forces.as_mut().map(|x| {
            fem_forces.extend_from_slice(x);
        });
        // M2
        //  - positioner
        fem_outputs
            .pop_these(vec![ios!(MCM2SmHexD)])
            .and_then(|mut hex_d| {
                hex_d.append(&mut vec![ios!(M2poscmd(vec![0f64; 42]))]);
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
        // WFSing
        if k >= wfs_delay {
            if let Some(ref mut atm) = gosm.atm {
                atm.secs = k as f64 / sampling_rate;
            }
            gosm.inputs(Some(vec![
                fem_outputs[ios!(OSSM1Lcl)].clone(),
                fem_outputs[ios!(MCM2Lcl6D)].clone(),
            ]))
            .unwrap()
            //gosm.inputs(fem_outputs.pop_these(ios!(OSSM1Lcl, MCM2Lcl6D))).unwrap()
            .step()
            .unwrap();
            if k % wfs_sample_rate == 0 {
                // Optical System
                //  - tip-tilt
                pzt_cmd = gosm
                    .outputs()
                    .and_then(|x| Into::<Option<Vec<f64>>>::into(x[0].clone()))
                    .map(|x| x.into_iter().map(|x| x as f32).collect::<Vec<f32>>())
                    .map(|y| wfs_2_m2rxy.solve(&mut y.into()).into())
                    .map(|m2_rxy: Vec<f32>| {
                        (rxy_2_stt
                            * na::DVector::from_iterator(14, m2_rxy.iter().map(|&x| x as f64)))
                        .as_slice()
                        .to_vec()
                    })
                    .map(|stt| os_tiptilt.in_step_out(Some(ios!(TTSP(vec![0f64]), TTFB(stt)))))
                    .unwrap()?;
            }
        }
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
        "{} sample wind loads played in [{}] ({}ms/step)",
        n_step,
        humantime::format_duration(eta),
        eta.as_millis() as f64 / n_step as f64
    );

    let (tilts, segments): (Vec<_>, Vec<_>) = data_log
        .into_iter()
        .collect::<Vec<f64>>()
        .chunks(23)
        .map(|data| {
            let (a, b) = data.split_at(2);
            (a.to_vec(), b.to_vec())
        })
        .unzip();
    let (segments_piston, _segments_tilts): (Vec<_>, Vec<_>) = segments
        .into_iter()
        .flatten()
        .collect::<Vec<f64>>()
        .chunks(21)
        .map(|data| {
            let (a, b) = data.split_at(7);
            (a.to_vec(), b.to_vec())
        })
        .unzip();

    Plot::from((
        tilts
            .into_iter()
            .flatten()
            .collect::<Vec<f64>>()
            .chunks(2)
            .enumerate()
            .map(|(k, tilts)| {
                (
                    k as f64 / sampling_rate,
                    tilts.iter().map(|x| x.to_mas()).collect::<Vec<f64>>(),
                )
            }),
        Some(
            Config::new()
                .filename("ns-opm-im_tilts.svg")
                .xaxis(Axis::new().label("Time [s]"))
                .yaxis(Axis::new().label("Tip-tilt [mas]")),
        ),
    ));

    Plot::from((
        segments_piston
            .into_iter()
            .flatten()
            .collect::<Vec<f64>>()
            .chunks(7)
            .enumerate()
            .map(|(k, pistons)| {
                (
                    k as f64 / sampling_rate,
                    pistons.iter().map(|x| x * 1e6).collect::<Vec<f64>>(),
                )
            }),
        Some(
            Config::new()
                .filename("ns-opm-im_pistons.svg")
                .xaxis(Axis::new().label("Time [s]"))
                .yaxis(Axis::new().label("Piston [micron]")),
        ),
    ));

    /*
    Plot::from((
        m1_axial
            .into_iter()
            .collect::<Vec<f64>>()
            .chunks(7)
            .enumerate()
            .map(|(k, y)| {
                (
                    k as f64 / sampling_rate,
                    y.iter().map(|y| y * 1e6).collect::<Vec<f64>>(),
                )
            }),
        Some(
            Config::new()
                .filename("ns-opm-im_m1-sfe-rms.svg")
                .xaxis(Axis::new().label("Time [s]"))
                .yaxis(Axis::new().label("M1 axial rms sfe [micron]")),
        ),
    ));*/

    /*
        let iter1: Vec<(f64, Vec<f64>)> = segments_tilts
            .into_iter()
            .flatten()
            .collect::<Vec<f64>>()
            .chunks(14)
            //.flat_map(|x| x.chunks(2).map(|x| x[0].hypot(x[1])).collect::<Vec<f64>>())
            .flat_map(|x| (0..7).map(|k| x[k].hypot(x[k + 7])).collect::<Vec<f64>>())
            .collect::<Vec<f64>>()
            .chunks(7)
            .enumerate()
            .map(|(k, tilts)| {
                (
                    k as f64 / sampling_rate,
                    tilts.iter().map(|x| x.to_mas()).collect::<Vec<f64>>(),
                )
            })
            .collect();
        let iter2: Vec<(f64, Vec<f64>)> = wfs_log
            .chunks(14)
            .flat_map(|x| x.chunks(2).map(|x| x[0].hypot(x[1])).collect::<Vec<f64>>())
            .collect::<Vec<f64>>()
            .chunks(7)
            .enumerate()
            .map(|(k, tilts)| {
                (
                    (k * wfs_sample_rate + wfs_delay) as f64 / sampling_rate,
                    tilts
                        .iter()
                        .map(|x| 0.25 * x.to_mas())
                        .collect::<Vec<f64>>(),
                )
            })
            .collect();

        let mut config = Config::new()
            .filename("ns-opm-im_m2-rxy.svg")
            .xaxis(Axis::new().label("Time [s]"))
            .yaxis(Axis::new().label("Segment tip-tilt mag. [mas]"));
        config.auto_range(vec![&iter1, &iter2]);
        <Combo as From<Complot>>::from((
            vec![Box::new(iter1.into_iter()), Box::new(iter2.into_iter())],
            vec![Kind::Plot, Kind::Scatter],
            Some(config),
        ));
    */
    Ok(())
}
