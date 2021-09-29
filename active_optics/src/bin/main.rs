use active_optics::QP;
use crseo::{
    calibrations,
    dos::GmtOpticalSensorModel,
    shackhartmann::{Diffractive, Geometric, WavefrontSensorBuilder},
    Builder, Calibration, ShackHartmann, GMT, SH48, SOURCE,
};
use dosio::{ios, Dos};
use skyangle::Conversion;
use std::{fs::File, time::Instant};

// Shack-Hartmann WFS kind
type WfsType = Diffractive;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // WFS Definition
    let n_aco_sensor = 3;
    let mut gosm_aco = GmtOpticalSensorModel::<ShackHartmann<WfsType>, SH48<WfsType>>::new(Some(
        SOURCE::new()
            .size(n_aco_sensor)
            .on_ring(6f32.from_arcmin())
            .fwhm(6f64),
    ))
    .gmt(GMT::new().m1("m1_eigen_modes", 329))
    .sensor(SH48::<WfsType>::new().n_sensor(n_aco_sensor))
    .build()?;
    println!("valid lenslets: {}", &gosm_aco.sensor.n_valid_lenslet());

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
    let file = File::open("../m1_eigen_modes.bin").unwrap();
    let data: Vec<(Vec<f64>, Vec<f64>)> = bincode::deserialize_from(file).unwrap();
    let (_, coefs2forces): (Vec<_>, Vec<_>) = data.into_iter().unzip();
    let n_actuators = vec![335, 335, 335, 335, 335, 335, 306];

    // Active optics control algorithm
    let mut aco = QP::<41, 41, 27, 271>::new(
        "SHAcO_qp_rhoP1e-3_kIp5.rs.pkl",
        (&dmat, 271),
        (&coefs2forces, n_actuators),
    )?
    .build();

    // Misalignment
    let mut m1_modes_buf = vec![vec![0f64; 329]; 7];
    m1_modes_buf[0][0] = 1e-6;
    m1_modes_buf[3][5] = 1e-6;
    m1_modes_buf[6][3] = -1e-6;
    let m1_modes = ios!(M1modes(
        m1_modes_buf.into_iter().flatten().collect::<Vec<f64>>()
    ));
    let mut m1_rbm_buf = vec![vec![0f64; 6]; 7];
    m1_rbm_buf[1][1] = 1e-6;
    let m1_rbm = ios!(OSSM1Lcl(
        m1_rbm_buf.into_iter().flatten().collect::<Vec<f64>>()
    ));
    let mut m2_rbm_buf = vec![0f64; 42];
    m2_rbm_buf[3] = 1e-6;
    let m2_rbm = ios!(MCM2Lcl6D(m2_rbm_buf));
    let y_valid = gosm_aco.in_step_out(Some(vec![m1_rbm, m2_rbm, m1_modes]))?;

    // reconstruction and integration
    let solution = aco.in_step_out(y_valid)?.unwrap();
    let ungain = -aco.controller_gain().recip();
    solution
        .into_iter()
        .inspect(|x| println!("{}", x.kind()))
        .filter_map(|x| -> Option<Vec<f64>> { x.into() })
        .for_each(|x| x.iter().for_each(|x| println!("{:+5.2}", ungain * x * 1e6)));

    Ok(())
}
