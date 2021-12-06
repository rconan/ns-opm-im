use active_optics::QP;
use complot::{Axis, Config, Plot};
use crseo::{
    calibrations, ceo,
    dos::GmtOpticalSensorModel,
    from_opticals,
    shackhartmann::{Diffractive, Geometric, WavefrontSensorBuilder},
    Builder, Calibration, OpticalSensitivities, ShackHartmann, GMT, SH48, SOURCE,
};
use dosio::{ios, Dos, IOTags, IOVec, IO};
use fem::{
    dos::{DiscreteModalSolver, DiscreteStateSpace, Exponential},
    FEM,
};
use geotrans::{Quaternion, Vector};
use indicatif::{ProgressBar, ProgressBarIter, ProgressIterator, ProgressStyle};
use m1_ctrl as m1;
use mount_ctrl as mount;
use nalgebra as na;
use serde_pickle as pkl;
use skyangle::Conversion;
use std::{fs::File, path::Path, time::Instant};
use windloading::WindLoads;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sim_duration = 15f64;
    let sampling_rate = 1e3;
    let wfs_delay = 10f64;
    let wfs_sample_delay = (sampling_rate * wfs_delay) as usize;
    let aco_delay = 6f64;
    let aco_sample_delay = (sampling_rate * aco_delay) as usize;
    let tiptilt_exposure_time = 5e-3;
    let tiptilt_sample_rate = (tiptilt_exposure_time * sampling_rate) as usize;
    let aco_exposure_time = 1f64;
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
        Path::new("data").join("b2019_0z_30az_os_7ms.wind_loads_1kHz_100-400s.pkl"),
    )?
    .range(0.0, sim_duration)
    .truss()?
    .build()?;
    println!("... in {}ms", now.elapsed().as_millis());

    // MOUNT
    let mut mnt_drives = mount::drives::Controller::new();
    let mut mnt_ctrl = mount::controller::Controller::new();
    // M1
    let mut m1_hardpoints = m1::hp_dynamics::Controller::new();
    let mut m1_load_cells = m1::hp_load_cells::Controller::new();
    let m1_actuator_sampling_rate = 100f64;
    let m1_actuators_update_rate = (sampling_rate / m1_actuator_sampling_rate) as usize;
    let mut m1_actuators = m1::actuators::M1ForceLoops::new();
    // M2
    let mut fsm_positionner = fsm::positionner::Controller::new();
    let mut fsm_piezostack = fsm::piezostack::Controller::new();
    // Optical System
    // - tipt-tilt
    let mut os_tiptilt = fsm::tiptilt::Controller::new();
    //FEM
    let sid = 1;
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
                .outputs(ios!(
                    M1Segment1AxialD,
                    M1Segment2AxialD,
                    M1Segment3AxialD,
                    M1Segment4AxialD,
                    M1Segment5AxialD,
                    M1Segment6AxialD,
                    M1Segment7AxialD
                ))
                .build()?,
        )
    };
    println!("... in {}ms", now.elapsed().as_millis());
    println!("{}", fem);
    /*{
        bincode::serialize_into(
            File::create("fem.bin")?,
            &(&m1_segments_surface_nodes, &fem),
        )?;
    }*/

    // CEO WFSS
    type WfsType = Diffractive;
    let n_tiptilt_sensor = 1;
    let mut gosm = GmtOpticalSensorModel::<ShackHartmann<WfsType>, SH48<WfsType>>::new(Some(
        SOURCE::new().size(n_tiptilt_sensor).fwhm(6f64),
    ))
    .gmt(GMT::new().m1("m1_eigen_modes", 329))
    .sensor(SH48::<WfsType>::new().n_sensor(n_tiptilt_sensor))
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
    // - M2 Rxy calibration matrix
    let mut gmt2wfs = Calibration::new(
        &gosm.gmt,
        &gosm.src,
        SH48::<Geometric>::new().n_sensor(n_tiptilt_sensor),
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
    // - poke matrix check
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

    // - M1 27 eigen modes calibration matrix
    let n_aco_sensor = 3;
    let mut gosm_aco = GmtOpticalSensorModel::<ShackHartmann<WfsType>, SH48<WfsType>>::new(Some(
        SOURCE::new()
            .size(n_aco_sensor)
            .on_ring(6f32.from_arcmin())
            .fwhm(6f64),
    ))
    .gmt(GMT::new().m1("m1_eigen_modes", 329))
    .sensor(SH48::<WfsType>::new().n_sensor(n_aco_sensor))
    /*        .atmosphere(crseo::ATMOSPHERE::new().ray_tracing(
        26.,
        520,
        0.,
        25.,
        Some("ns-opm-im_atm.bin".to_string()),
        Some(8),
    ))*/
    .build()?;

    // WFS Calibration
    use calibrations::Mirror::*;
    use calibrations::Segment::*;
    let mut segments = vec![vec![Txyz(1e-6, None), Rxyz(1e-6, None)]; 6];
    segments.append(&mut vec![vec![Txyz(1e-6, None), Rxyz(1e-6, Some(0..2))]]);
    let dmat: Vec<_> = vec![
        (vec![M1], segments.clone()),
        (vec![M2], segments),
        (vec![M1MODES], vec![vec![Modes(1e-6, 0..27)]; 7]),
    ]
    .into_iter()
    .map(|(mirror, segments)| -> Vec<f32> {
        let gmt_calib = GMT::new().m1("m1_eigen_modes", 27).build().unwrap();
        let wfs_calib = SH48::<Geometric>::new().n_sensor(n_aco_sensor);
        let src_calib = wfs_calib
            .guide_stars(Some(
                SOURCE::new().size(n_aco_sensor).on_ring(6f32.from_arcmin()),
            ))
            .build()
            .unwrap();
        // WFS Calibration
        println!("M1 RBM calibration");
        let mut calib = Calibration::new(&gmt_calib, &src_calib, wfs_calib.clone());
        let now = Instant::now();
        calib.calibrate(
            mirror,
            segments,
            calibrations::ValidLensletCriteria::OtherSensor(&mut gosm_aco.sensor),
        );
        println!(
            "GMT WFS calibration [{}x{}] in {}s",
            &calib.n_data,
            &calib.n_mode,
            now.elapsed().as_secs()
        );
        calib.poke.into()
    })
    .flatten()
    .map(|x| x as f64)
    .collect();

    // Eigen modes and eigen coefs to forces matrix
    let file = File::open("m1_eigen_modes.bin").unwrap();
    let data: Vec<(Vec<f64>, Vec<f64>)> = bincode::deserialize_from(file).unwrap();
    let (_, coefs2forces): (Vec<_>, Vec<_>) = data.into_iter().unzip();
    let n_actuators = vec![335, 335, 335, 335, 335, 335, 306];

    // Active optics control algorithm
    let mut aco = QP::<41, 41, 27, 271>::new(
        "active_optics/SHAcO_qp_rhoP1e-3_kIp5.rs.pkl",
        (&dmat, 271),
        (&coefs2forces, n_actuators),
    )?
    //    .convergence_tolerances((1e-6, 1e-4))
    //    .verbose()
    .as_m1_actuator_forces()
    .build();

    // Eigen modes and eigen coefs to forces matrix
    let file = File::open("m1_eigen_modes.bin").unwrap();
    let eigens_data: Vec<(Vec<f64>, Vec<f64>)> = bincode::deserialize_from(file).unwrap();
    let (eigens, coefs2forces) = &eigens_data[sid - 1];
    let nodes = &m1_segments_surface_nodes[0];
    let n_node = nodes.len() / 3;
    let n_eigen_mode = eigens.len() / n_node;
    let m1_s1_eigens = na::DMatrix::from_column_slice(n_node, n_eigen_mode, eigens);

    let optical_sensitivities = OpticalSensitivities::load()?;
    let rxy_2_stt = &optical_sensitivities[3].m2_rxy()?;
    let opticals = from_opticals(&optical_sensitivities[1..4]);

    // I/O initialization
    let m1_s1_coefs = na::DVector::from_column_slice(&[5e-6, -2.5e-6, 3e-6]);
    let m1_s1_forces = na::DMatrix::from_column_slice(
        coefs2forces.len() / n_eigen_mode,
        n_eigen_mode,
        &coefs2forces,
    )
    .columns(0, m1_s1_coefs.nrows())
        * &m1_s1_coefs;
    dbg!(m1_s1_forces.shape());
    /*
        let mut m1_actuators_forces = Some(ios!(
            M1ActuatorsSegment1(m1_s1_forces.as_slice().to_vec()),
            M1ActuatorsSegment2(vec![0f64; 335]),
            M1ActuatorsSegment3(vec![0f64; 335]),
            M1ActuatorsSegment4(vec![0f64; 335]),
            M1ActuatorsSegment5(vec![0f64; 335]),
            M1ActuatorsSegment6(vec![0f64; 335]),
            M1ActuatorsSegment7(vec![0f64; 306])
        ));
    */
    let mut m1_bending_modes = ios!(
        M1S1BMcmd(vec![0f64; 335]),
        //M1S1BMcmd(m1_s1_forces.as_slice().to_vec()),
        M1S2BMcmd(vec![0f64; 335]),
        M1S3BMcmd(vec![0f64; 335]),
        M1S4BMcmd(vec![0f64; 335]),
        M1S5BMcmd(vec![0f64; 335]),
        M1S6BMcmd(vec![0f64; 335]),
        M1S7BMcmd(vec![0f64; 306])
    );

    let mut m1_hardpoints_forces = Some(vec![ios!(OSSHarpointDeltaF(vec![0f64; 42]))]);
    let mut m1_actuators_forces = None;
    let mut pzt_cmd = Some(vec![ios!(TTcmd(vec![0f64; 21]))]);

    let mut fem_outputs = fem.zeroed_outputs();

    let norm =
        |x: &[f64]| -> f64 { (x.iter().map(|x| x * x).sum::<f64>() / x.len() as f64).sqrt() };
    dbg!(norm(m1_s1_forces.as_slice()));

    let n_step = wind_loading.n_sample;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m2_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut data_log = Vec::<f64>::with_capacity(n_step * 23);

    let pb = ProgressBar::new(n_step as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{duration_precise}] {bar:60.cyan/blue} {elapsed_precise} {msg}")
            .progress_chars("=|-"),
    );

    let now = Instant::now();
    (0..n_step).progress_with(pb).for_each(|k| {
        let mut fem_forces: Vec<IO<Vec<f64>>> = wind_loading.outputs().unwrap_or_default();
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
        // - actuators
        if k % m1_sample_rate == 0 {
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
        if k >= wfs_sample_delay {
            let bm_coefs: Vec<f64> = ios!(
                M1Segment1AxialD,
                M1Segment2AxialD,
                M1Segment3AxialD,
                M1Segment4AxialD,
                M1Segment5AxialD,
                M1Segment6AxialD,
                M1Segment7AxialD
            )
            .into_iter()
            .filter_map(|io| Option::<Vec<f64>>::from(&fem_outputs[io]))
            .zip(eigens_data.iter())
            .map(|(m1_segment_figure, (eigens, _))| {
                {
                    na::DMatrix::from_column_slice(
                        m1_segment_figure.len(),
                        eigens.len() / m1_segment_figure.len(),
                        eigens,
                    )
                    //                    .columns(0, 3)
                    .transpose()
                        * na::DVector::from_column_slice(&m1_segment_figure)
                }
                .as_slice()
                .to_vec()
            })
            .flat_map(|mut x| {
                if x.len() < 329 {
                    x.append(&mut vec![0f64; 329 - x.len()])
                }
                x
            })
            .collect();
            if let Some(ref mut atm) = gosm.atm {
                atm.secs = k as f64 / sampling_rate;
            }
            let os_gmt_state = Some(vec![
                fem_outputs[ios!(OSSM1Lcl)].clone(),
                fem_outputs[ios!(MCM2Lcl6D)].clone(),
                ios!(M1modes(bm_coefs.clone())),
            ]);
            gosm.inputs(os_gmt_state.clone())
                .unwrap()
                //gosm.inputs(fem_outputs.pop_these(ios!(OSSM1Lcl, MCM2Lcl6D))).unwrap()
                .step()
                .unwrap();
            if k % tiptilt_sample_rate == 0 {
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
                    .unwrap()
                    .unwrap();
            }
            if k >= aco_sample_delay + wfs_sample_delay {
                gosm_aco.inputs(os_gmt_state).unwrap().step().unwrap();
                if k % aco_sample_rate == 0 {
                    // getting the slopes
                    /*let y_valid = gosm_aco.outputs().unwrap();
                    let slopes: Vec<f64> = y_valid[ios!(SensorData)]
                        .as_ref()
                        .and_then(|x| x.into())
                        .unwrap()
                        .to_vec();*/
                    // reconstruction and integration
                    m1_bending_modes = aco.in_step_out(gosm_aco.outputs()).unwrap().unwrap();
                    /*
                                            let sol_m1_s1_forces: &[f64] = solution[ios!(M1S1BMcmd)]
                                                .as_ref()
                                                .and_then(|x| x.into())
                                                .unwrap();
                                            // Set point forces on M1 S1 F0
                                            // Aco control provide ~ -k * F0
                                            // the command applied is F0 - (~ k *F0)
                                            let forces: Vec<f64> = m1_s1_forces
                                                .as_slice()
                                                .iter()
                                                .zip(sol_m1_s1_forces)
                                                .map(|(a, b)| *a + *b)
                                                .collect();
                                            /*                    println!(
                                                                    r##"
                                            {:-<12}
                                            Iteration step #{:05}:
                                             - WFS norm: {:.3}
                                             - BM coefs: {:.3?}
                                             - M1 actuators forces: {:.6}
                                            {:-<12}
                                            "##,
                                                                    "",
                                                                    k,
                                                                    norm(&slopes) * 1e6,
                                                                    bm_coefs
                                                                        .iter()
                                                                        .take(3)
                                                                        .map(|x| x * 1e6)
                                                                        .collect::<Vec<f64>>(),
                                                                    norm(&forces),
                                                                    ""
                                                                );*/
                                            m1_bending_modes[0] = ios!(M1S1BMcmd(forces));
                    */
                }
            }
        };
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
    });
    /*
            .filter_map(|m1_s1_coefs_from_figure: Option<Vec<f64>>| m1_s1_coefs_from_figure)
            .last()
            .map(|m1_s1_coefs_from_figure| {
                dbg!(m1_s1_coefs_from_figure
                    .as_slice()
                    .iter()
                    .take(3)
                    .map(|x| x * 1e6)
                    .collect::<Vec<f64>>());
            });
    */

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

    Ok(())
}
