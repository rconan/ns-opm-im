use crseo::{
    calibrations,
    dos::GmtOpticalSensorModel,
    from_opticals,
    imaging::NoiseDataSheet,
    sensitivities::Bin,
    shackhartmann::{Diffractive, Geometric, WavefrontSensorBuilder},
    Builder, Calibration, OpticalSensitivities, ShackHartmann, SH48,
};
use dosio::{ios, Dos, IOVec, IO};
use fem::{dos::DiscreteStateSpace, FEM};
use indicatif::{ProgressBar, ProgressStyle};
use m1_ctrl as m1;
use mount_ctrl as mount;
use nalgebra as na;
use ns_opm_im::complot::{Axis, Combo, Complot, Config, Kind, Plot};
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
    .range(0.0, 200.0)
    .truss()?
    .build()?;

    // MOUNT CONTROL
    let mut mnt_drives = mount::drives::Controller::new();
    let mut mnt_ctrl = mount::controller::Controller::new();

    // M1
    let mut m1_hardpoints = m1::hp_load_cells::Controller::new();
    let mut m1_force_loops = m1::actuators::M1ForceLoops::new();

    // FEM
    let m1_axial_ios = ios!(
        M1Segment1AxialD,
        M1Segment2AxialD,
        M1Segment3AxialD,
        M1Segment4AxialD,
        M1Segment5AxialD,
        M1Segment6AxialD,
        M1Segment7AxialD
    );
    let sampling_rate = 1e3;
    let m1_rbm_io = ios!(OSSM1Lcl);
    let m2_rbm_io = ios!(MCM2Lcl6D);
    let mut fem = {
        let fem = FEM::from_env()?;
        DiscreteStateSpace::from(fem)
    }
    .sampling(sampling_rate)
    .proportional_damping(2. / 100.)
    .inputs_from(&wind_loading)
    .inputs_from(&mnt_drives)
    //    .inputs(vec![jar::OSSM1Lcl6F::io()])
    .outputs(vec![m1_rbm_io.clone(), m2_rbm_io.clone()])
    .outputs(ios!(
        OSSAzEncoderAngle,
        OSSElEncoderAngle,
        OSSRotEncoderAngle
    ))
    .inputs(ios!(
        M1ActuatorsSegment1,
        M1ActuatorsSegment2,
        M1ActuatorsSegment3,
        M1ActuatorsSegment4,
        M1ActuatorsSegment5,
        M1ActuatorsSegment6,
        M1ActuatorsSegment7
    ))
    .outputs(vec![ios!(OSSHardpointD)])
    .outputs(m1_axial_ios.clone())
    .build()?;
    println!(
        "Dynamic discrete FEM: {} 2x2 state space model",
        fem.state_space.len()
    );

    let mut mount_drives_forces = Some(ios!(
        OSSAzDriveTorque(vec![0f64; 12]),
        OSSElDriveTorque(vec![0f64; 4]),
        OSSRotDriveTorque(vec![0f64; 4])
    ));

    // OPTICAL SENSITIVITIES
    /*    let opticals = from_opticals(
        &{
            if let Ok(senses) = OpticalSensitivities::load() {
                senses
            } else {
                OpticalSensitivities::compute(None)?.dump()?
            }
        }[1..4],
    );*/
    let opticals = from_opticals(&{ OpticalSensitivities::load().unwrap() }[1..4]);
    println!("opticals: {:?}", opticals.shape());

    //    let data_rm: |…| -> {unknown} = |data: &mut Vec<IO<Vec<f64>>>, p: IO<()>| {data.remove(data.iter().position(|x: &{unknown}| *x == p).unwrap());};

    // CEO WFSS
    let exposure_time = 10.;
    let n_sensor = 1;
    let wfs_sample_rate = (exposure_time * sampling_rate) as usize;
    let wfs_delay = sampling_rate as usize * 10;
    type WFS_TYPE = Diffractive;
    let mut gosm = GmtOpticalSensorModel::<ShackHartmann<WFS_TYPE>, SH48<WFS_TYPE>>::new()
        .sensor(SH48::<WFS_TYPE>::new().n_sensor(n_sensor))
        .atmosphere(crseo::ATMOSPHERE::new().ray_tracing(
            26.,
            520,
            0.,
            25.,
            Some("ns-opm-im_atm.bin".to_string()),
            Some(8),
        ))
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
        .in_step_out(Some(vec![m2_rbm.clone()]))
        .map(|x| Into::<Option<Vec<f64>>>::into(x.unwrap()[0].clone()))
        .unwrap()
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

    let mut data_log = Vec::<f64>::with_capacity(wind_loading.n_sample * 23);
    let mut m2_log = Vec::<f64>::with_capacity(wind_loading.n_sample * 14);
    let mut wfs_log = Vec::<f64>::with_capacity(wind_loading.n_sample * 14 / wfs_sample_rate);
    //    let mut mount_log = Vec::<f64>::with_capacity(wind_loading.n_sample);
    let mut m1_axial = Vec::<f64>::with_capacity(wind_loading.n_sample * 7);

    let pb = ProgressBar::new(wind_loading.n_sample as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:60.cyan/blue} {pos:>7}/{len:7}")
            .progress_chars("=|-"),
    );
    let now = Instant::now();
    let mut m1_actuators: Option<Vec<IO<Vec<f64>>>> = None;
    let mut k = 0;
    while let Some(mut fem_forces) = wind_loading.outputs() {
        pb.inc(1);
        // FEM
        if let Some(x) = mount_drives_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        if let Some(ref mut value) = m1_actuators {
            fem_forces.append(&mut value.clone());
        } else {
            fem_forces.append(&mut ios!(
                M1ActuatorsSegment1(vec![0f64; 335]),
                M1ActuatorsSegment2(vec![0f64; 335]),
                M1ActuatorsSegment3(vec![0f64; 335]),
                M1ActuatorsSegment4(vec![0f64; 335]),
                M1ActuatorsSegment5(vec![0f64; 335]),
                M1ActuatorsSegment6(vec![0f64; 335]),
                M1ActuatorsSegment7(vec![0f64; 306])
            ));
        };
        /*Option::<Vec<f64>>::from(&fem_forces[ios!(M1ActuatorsSegment1)])
            .map(|x| x.iter().map(|x| x * x).sum::<f64>().sqrt())
        .map(|x| mount_log.push(x));*/
        let mut fem_outputs = fem
            .in_step_out(Some(fem_forces))?
            .ok_or("FEM output is empty")?;
        m1_axial.extend(
            m1_axial_ios
                .iter()
                .map(|io| fem_outputs[io].mean())
                .collect::<Vec<f64>>()
                .into_iter(),
        );
        // WFSing
        let m2_rxy: Option<Vec<f32>> = if k >= wfs_delay {
            if let Some(ref mut atm) = gosm.atm {
                atm.secs = k as f64 / sampling_rate;
            }
            gosm.inputs(Some(vec![
                fem_outputs[ios!(OSSM1Lcl)].clone(),
                fem_outputs[ios!(MCM2Lcl6D)].clone(),
            ]))?
            //gosm.inputs(fem_outputs.pop_these(ios!(OSSM1Lcl, MCM2Lcl6D)))?
            .step()?;
            if k % wfs_sample_rate == 0 {
                let y = gosm
                    .outputs()
                    .map(|x| Into::<Option<Vec<f64>>>::into(x[0].clone()))
                    .unwrap()
                    .map(|x| x.into_iter().map(|x| x as f32).collect::<Vec<f32>>())
                    .unwrap();
                Some(wfs_2_m2rxy.solve(&mut y.into()).into())
            } else {
                None
            }
        } else {
            None
        };
        // MOUNT CONTROLLER & DRIVES
        let mount_encoders = vec![
            fem_outputs[ios!(OSSElEncoderAngle)].clone(),
            fem_outputs[ios!(OSSAzEncoderAngle)].clone(),
            fem_outputs[ios!(OSSRotEncoderAngle)].clone(),
        ];
        mount_drives_forces = mnt_ctrl
            .in_step_out(Some(mount_encoders.clone()))?
            .map(|mut x| {
                x.extend_from_slice(&mount_encoders);
                mnt_drives.in_step_out(Some(x))
            })
            .unwrap()?;
        // M1 HARDPOINTS/ACTUATORS FORCE LOOP
        if k % 10 == 0 {
            let m1_hp = vec![
                <Vec<IO<Vec<f64>>> as IOVec>::pop_this(&mut fem_outputs, ios!(OSSHardpointD))
                    .unwrap(),
                ios!(M1HPCmd(vec![0f64; 42])),
            ];
            let mut hps = m1_hardpoints.in_step_out(Some(m1_hp))?;
            if let Some(ref mut hps) = hps {
                hps.append(&mut vec![
                    ios!(M1S1BMcmd(vec![0f64; 27])),
                    ios!(M1S2BMcmd(vec![0f64; 27])),
                    ios!(M1S3BMcmd(vec![0f64; 27])),
                    ios!(M1S4BMcmd(vec![0f64; 27])),
                    ios!(M1S5BMcmd(vec![0f64; 27])),
                    ios!(M1S6BMcmd(vec![0f64; 27])),
                    ios!(M1S7BMcmd(vec![0f64; 27])),
                ]);
            };
            m1_actuators = m1_force_loops.in_step_out(hps)?;
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

        // M2 Rxy logging
        m2_log.append(
            &mut Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)])
                .map(|x| {
                    x.chunks(6)
                        .map(|x| x[3..5].to_vec())
                        .flatten()
                        .collect::<Vec<f64>>()
                })
                .unwrap(),
        );
        // WFS M2 Rxy estimate logging
        if let Some(values) = m2_rxy {
            wfs_log.extend(values.iter().map(|&x| x as f64))
        }
    }
    pb.finish();
    let eta = now.elapsed();
    println!(
        "{} sample wind loads played in [{}] ({}ms/step)",
        wind_loading.n_sample,
        humantime::format_duration(eta),
        eta.as_millis() as f64 / wind_loading.n_sample as f64
    );

    /*
        println!("m2_log: {}", m2_log.len());
        println!("wfs_log: {:#?}", wfs_log);
        let u_f64 = m2_log.chunks(14).last().unwrap().to_vec();
        let u_f32 = u_f64.iter().map(|&x| x as f32).collect::<Vec<f32>>();
        println!("M2 rxy");
        u_f64
            .chunks(2)
            .enumerate()
            .for_each(|x| println!("#{}: [{:+0.1},{:+0.1}]", 1 + x.0, x.1[0], x.1[1]));
        gosm.gmt.reset();
        let y = gosm
            .in_step_out(Some(vec![ios!(MCM2Lcl6D(u_f64))]))
            .map(|x| Into::<Option<Vec<f64>>>::into(x.unwrap()[0].clone()))
            .unwrap()
            .map(|x| x.into_iter().map(|x| x as f32).collect::<Vec<f32>>())
            .unwrap();
        let a: Vec<f32> = wfs_2_m2rxy.solve(&mut y.into()).into();
        a.into_iter()
            .map(|x| x * 1e6)
            .collect::<Vec<f32>>()
            .chunks(2)
            .enumerate()
            .for_each(|x| println!("#{}: [{:+0.1},{:+0.1}]", 1 + x.0, x.1[0], x.1[1]));
        //    println!("WFS M2 rxy: {:#?}", a);
    */
    let (tilts, segments): (Vec<_>, Vec<_>) = data_log
        .into_iter()
        .collect::<Vec<f64>>()
        .chunks(23)
        .map(|data| {
            let (a, b) = data.split_at(2);
            (a.to_vec(), b.to_vec())
        })
        .unzip();
    let (segments_piston, segments_tilts): (Vec<_>, Vec<_>) = segments
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
                .yaxis(Axis::new().label("Piston [nm]")),
        ),
    ));

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
    ));

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

    Ok(())
}
