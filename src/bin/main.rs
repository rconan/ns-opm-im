use active_optics::QP;
use complot::{Axis, Config, Plot};
use crseo::{
    calibrations,
    dos::GmtOpticalSensorModel,
    from_opticals,
    shackhartmann::{Diffractive, Geometric, WavefrontSensorBuilder},
    Builder, Calibration, OpticalSensitivities, ShackHartmann, GMT, SH48, SOURCE,
};
use dosio::{ios, Dos, IOVec, IO};
use fem::{
    dos::{DiscreteModalSolver, DiscreteStateSpace, Exponential},
    FEM,
};
use geotrans::{Quaternion, Vector};
use indicatif::{ProgressBar, ProgressStyle};
use m1_ctrl as m1;
use mount_ctrl as mount;
use nalgebra as na;
use skyangle::Conversion;
use std::{fs::File, io::BufWriter, path::Path, time::Instant};
use windloading::WindLoads;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sim_duration = 60f64;
    let sampling_rate = 1e3;
    let wfs_delay = 10f64;
    let wfs_sample_delay = (sampling_rate * wfs_delay) as usize;
    let aco_delay = 1f64;
    let aco_sample_delay = (sampling_rate * aco_delay) as usize;
    let tiptilt_exposure_time = 5e-3;
    let tiptilt_sample_rate = (tiptilt_exposure_time * sampling_rate) as usize;
    let aco_exposure_time = 5f64;
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
    .m1_segments()?
    .m1_cell()?
    .m2_asm_topend()?
    .m2_segments_into(ios!(MCM2LclForce6F))?
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
    // - tipt-tilt
    let mut os_tiptilt = fsm::tiptilt::Controller::new();
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
                .truncate_hankel_singular_values(1e-5)
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

    // CEO WFSS
    let n_sensor = 1;
    type WfsType = Diffractive;
    let mut gosm = GmtOpticalSensorModel::<ShackHartmann<WfsType>, SH48<WfsType>>::new(Some(
        SOURCE::new().size(n_sensor).fwhm(6f64),
    ))
    .gmt(GMT::new().m1("m1_eigen_modes", 329))
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
    let wfs_2_m2rxy = gmt2wfs.qr();

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
    let eigens_data: Vec<(Vec<f64>, Vec<f64>)> = bincode::deserialize_from(file).unwrap();
    let (eigens, coefs2forces): (Vec<_>, Vec<_>) = eigens_data.into_iter().unzip();
    let n_actuators = vec![335, 335, 335, 335, 335, 335, 306];
    let modes_t: Vec<_> = eigens
        .iter()
        .zip(n_actuators.iter().map(|&n| n - 6))
        .map(|(x, n)| na::DMatrix::from_column_slice(x.len() / n, n, x).transpose())
        .collect();

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
    let mut m1_bending_modes = ios!(
        M1S1BMcmd(vec![0f64; 335]),
        M1S2BMcmd(vec![0f64; 335]),
        M1S3BMcmd(vec![0f64; 335]),
        M1S4BMcmd(vec![0f64; 335]),
        M1S5BMcmd(vec![0f64; 335]),
        M1S6BMcmd(vec![0f64; 335]),
        M1S7BMcmd(vec![0f64; 306])
    );
    let mut m1_hardpoints_forces = Some(vec![ios!(OSSHarpointDeltaF(vec![0f64; 42]))]);
    let mut m1_rbm_cmd: Option<Vec<IO<Vec<f64>>>> = None;
    let mut m2_pos_cmd = vec![ios!(M2poscmd(vec![0f64; 42]))];
    let mut pzt_cmd = Some(vec![ios!(TTcmd(vec![0f64; 21]))]);
    let mut fem_outputs = fem.zeroed_outputs();

    let n_step = wind_loading.n_sample;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m1_bm_logs = Vec::<Vec<f64>>::with_capacity(n_step * (329 * 6 + 300));
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
        .filter_map(|io| fem_outputs.pop_this(io))
        .filter_map(|x| Option::<Vec<f64>>::from(x))
        .zip(&m1_segments_surface_nodes)
        .zip(m1_rbm.as_ref().unwrap().chunks(6))
        .map(|((shape, nodes), rbm)| {
            let (t_xyz, r_xyz) = rbm.split_at(3);
            let rxyz_surface = {
                let q = Quaternion::unit(r_xyz[2], Vector::k())
                    * Quaternion::unit(r_xyz[1], Vector::j())
                    * Quaternion::unit(r_xyz[0], Vector::i());
                let trans_nodes: Vec<f64> = nodes
                    .chunks(3)
                    .flat_map(|x| {
                        let p: Quaternion = Vector::from(x).into();
                        let w: Quaternion = &q * p * &q.complex_conjugate();
                        Vec::<f64>::from(Vector::from(w.vector_as_slice()))
                    })
                    .collect();
                trans_nodes
                    .chunks(3)
                    .map(|x| x[2])
                    .zip(nodes.chunks(3).map(|x| x[2]))
                    .map(|(z_rbm, z)| z_rbm - z)
                    .collect::<Vec<f64>>()
            };
            shape
                .iter()
                .zip(rxyz_surface.iter())
                .map(|(a, r)| a - r - t_xyz[2])
                .collect::<Vec<f64>>()
        })
        .zip(&modes_t)
        .map(|(m1_segment_figure, eigens)| {
            { eigens * na::DVector::from_column_slice(&m1_segment_figure) }
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
        //        bm_coefs[0] += 1e-6;
        m1_bm_logs.push(bm_coefs.clone());
        // WFSing
        if k >= wfs_sample_delay {
            let os_gmt_state = Some(vec![
                fem_outputs[ios!(OSSM1Lcl)].clone(),
                fem_outputs[ios!(MCM2Lcl6D)].clone(),
                ios!(M1modes(bm_coefs)),
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
                    .unwrap()?;
            }
            let l = k - wfs_sample_delay;
            if l >= aco_sample_delay {
                gosm_aco.inputs(os_gmt_state).unwrap().step().unwrap();
                if l % aco_sample_rate == 0 {
                    m1_bending_modes = aco.in_step_out(gosm_aco.outputs()).unwrap().unwrap();
                    m2_pos_cmd = vec![m1_bending_modes.pop().unwrap()];
                    m1_rbm_cmd = m1_bending_modes.pop().map(|x| vec![x]);
                }
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

        if k % 1000 == 0 {
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
            let mut wfe_rms = vec![];
            for mut bm in m1_bm_logs.clone() {
                let mut b_rss = vec![];
                for n in [329, 329, 329, 329, 329, 329, 300] {
                    let b = bm.drain(..n);
                    b_rss.push(b.into_iter().map(|x| x * x).sum::<f64>().sqrt());
                }
                wfe_rms.push(b_rss);
            }

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

            Plot::from((
                wfe_rms.into_iter().enumerate().map(|(k, wfe_rms)| {
                    (
                        k as f64 / sampling_rate,
                        wfe_rms.iter().map(|x| x * 1e6).collect::<Vec<f64>>(),
                    )
                }),
                Some(
                    Config::new()
                        .filename("ns-opm-im_wfes.svg")
                        .xaxis(Axis::new().label("Time [s]"))
                        .yaxis(Axis::new().label("WFE RMS [micron]")),
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
                            .filename(&format!("ns-opm-im_m1-s{}-rbm.svg", sid + 1))
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
                            .filename(&format!("ns-opm-im_m2-s{}-rbm.svg", sid + 1))
                            .xaxis(Axis::new().label("Time [s]"))
                            .yaxis(Axis::new().label("WFE RMS [micron]")),
                    ),
                ));
            });
        }
    }
    pb.finish();
    let eta = now.elapsed();
    println!(
        "{} sample wind loads played in [{}] ({:.3}ms/step)",
        n_step,
        humantime::format_duration(eta),
        eta.as_millis() as f64 / n_step as f64
    );

    let file = File::create("ns-opm-im_gmt-logs.pkl")?;
    let mut buf = BufWriter::with_capacity(1_000_000, file);
    serde_pickle::to_writer(&mut buf, &(m1_logs, m2_logs, m1_bm_logs), true)?;

    Ok(())
}
