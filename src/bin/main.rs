use crseo::{dos::GmtOpticalSensorModel, from_opticals, sensitivities::Bin, OpticalSensitivities};
use crseo::{
    imaging::NoiseDataSheet, shackhartmann::Diffractive as WFS_TYPE,
    shackhartmann::WavefrontSensorBuilder, Builder, ShackHartmann, SH48,
};
use dosio::{io::jar, ios, Dos, IOVec, IO};
use fem::{dos::DiscreteStateSpace, FEM};
use indicatif::ProgressBar;
use m1_ctrl as m1;
use mount_ctrl as mount;
use nalgebra as na;
use ns_opm_im::complot::{Axis, Config, Plot};
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
    .range(0.0, 20.0)
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
    let m1_rbm_io = jar::OSSM1Lcl::io();
    let m2_rbm_io = jar::MCM2Lcl6D::io();
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

    let mut data_log = Vec::<f64>::with_capacity(wind_loading.n_sample * 23);
    //    let mut mount_log = Vec::<f64>::with_capacity(wind_loading.n_sample);
    let mut m1_axial = Vec::<f64>::with_capacity(wind_loading.n_sample * 7);

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

    let exposure_time = 1.0;
    let wfs_sample_rate = (exposure_time * sampling_rate) as usize;
    let mut gosm = GmtOpticalSensorModel::<ShackHartmann<WFS_TYPE>, SH48<WFS_TYPE>>::new()
        .sensor(
            SH48::<WFS_TYPE>::new()
                .n_sensor(1)
                .detector_noise_specs(NoiseDataSheet::new(exposure_time)),
        )
        .build()?;
    println!("M1 mode: {}", gosm.gmt.get_m1_mode_type());
    println!("M2 mode: {}", gosm.gmt.get_m2_mode_type());
    println!("GS band: {}", gosm.src.get_photometric_band());

    let pb = ProgressBar::new(wind_loading.n_sample as u64);
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
        gosm.inputs(Some(vec![
            fem_outputs[ios!(OSSM1Lcl)].clone(),
            fem_outputs[ios!(MCM2Lcl6D)].clone(),
        ]))?
        //gosm.inputs(fem_outputs.pop_these(ios!(OSSM1Lcl, MCM2Lcl6D)))?
        .step()?;
        let sensor_data = if k > 0 && k % wfs_sample_rate == 0 {
            gosm.outputs()
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
    }
    pb.finish();
    let eta = now.elapsed().as_millis();
    println!(
        "{} sample wind loads played in {}ms ({}ms/step)",
        wind_loading.n_sample,
        eta,
        eta as f64 / wind_loading.n_sample as f64
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

    Ok(())
}
