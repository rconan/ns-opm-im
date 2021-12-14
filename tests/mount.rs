use dosio::{ios, Dos};
use fem::{
    dos::{DiscreteStateSpace, Exponential},
    FEM,
};
use mount_ctrl as mount;
use serde_pickle as pkl;
use skyangle::Conversion;
use std::fs::File;

#[test]
fn mount_zeros() {
    // MOUNT CONTROL
    let mut mnt_drives = mount::drives::Controller::new();
    let mut mnt_ctrl = mount::controller::Controller::new();
    //FEM
    let sampling_rate = 1e3;
    let mut fem = {
        let fem = FEM::from_env().unwrap();
        println!("FEM eigen frequencies: {:?}", &fem.eigen_frequencies[..5]);
        DiscreteStateSpace::from(fem)
    }
    .sampling(sampling_rate)
    .proportional_damping(2. / 100.)
    .inputs_from(&[&mnt_drives])
    .outputs(ios!(OSSM1Lcl, MCM2Lcl6D))
    .outputs(ios!(
        OSSAzEncoderAngle,
        OSSElEncoderAngle,
        OSSRotEncoderAngle
    ))
    .build::<Exponential>()
    .unwrap();

    let mut mount_drives_forces = Some(ios!(
        OSSAzDriveTorque(vec![0f64; 12]),
        OSSElDriveTorque(vec![0f64; 4]),
        OSSRotDriveTorque(vec![0f64; 4])
    ));

    for _ in 0..3 {
        let fem_outputs = fem.in_step_out(mount_drives_forces).unwrap().unwrap();

        let mount_encoders = vec![
            fem_outputs[ios!(OSSElEncoderAngle)].clone(),
            fem_outputs[ios!(OSSAzEncoderAngle)].clone(),
            fem_outputs[ios!(OSSRotEncoderAngle)].clone(),
        ];
        mount_drives_forces = mnt_ctrl
            .in_step_out(Some(mount_encoders.clone()))
            .unwrap()
            .map(|mut x| {
                x.extend_from_slice(&mount_encoders);
                mnt_drives.in_step_out(Some(x))
            })
            .unwrap()
            .unwrap();
        let el = Option::<Vec<f64>>::from(
            &mount_drives_forces.as_mut().unwrap()[ios!(OSSElDriveTorque)],
        )
        .unwrap();
        println!(
            "Elevation: {}",
            el.iter().map(|x| x * x).sum::<f64>().sqrt()
        );
    }

    let fem_outputs = fem.outputs().unwrap();
    let m1_rbm = Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)]).unwrap();
    let m2_rbm = Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).unwrap();
    m1_rbm.chunks(7).enumerate().for_each(|(sid, rbm)| {
        let (t_xyz, r_xyz) = rbm.split_at(3);
        println!("M1 Segment #{}:", sid + 1);
        let tn = 1e9 * t_xyz.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(tn < 1e-1);
        println!(" - T: {:9.6}nm", tn);
        let rn = r_xyz.iter().map(|x| x * x).sum::<f64>().sqrt().to_mas();
        assert!(rn < 1e-1);
        println!(" - R: {:9.6}mas", rn);
    });
    m2_rbm.chunks(7).enumerate().for_each(|(sid, rbm)| {
        let (t_xyz, r_xyz) = rbm.split_at(3);
        println!("M2 Segment #{}:", sid + 1);
        let tn = 1e9 * t_xyz.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(tn < 1e-1);
        println!(
            " - T: {:9.6}nm",
            1e9 * t_xyz.iter().map(|x| x * x).sum::<f64>().sqrt()
        );
        let rn = r_xyz.iter().map(|x| x * x).sum::<f64>().sqrt().to_mas();
        assert!(rn < 1e-1);
        println!(
            " - R: {:9.6}mas",
            r_xyz.iter().map(|x| x * x).sum::<f64>().sqrt().to_mas()
        );
    });
}

#[test]
fn mount_constants() {
    // MOUNT CONTROL
    let mut mnt_drives = mount::drives::Controller::new();
    let mut mnt_ctrl = mount::controller::Controller::new();
    //FEM
    let sampling_rate = 1e3;
    let mut fem = {
        let fem = FEM::from_env().unwrap();
        println!("FEM\n{}", fem);
        DiscreteStateSpace::from(fem)
    }
    .sampling(sampling_rate)
    .proportional_damping(2. / 100.)
    .inputs_from(&[&mnt_drives])
    .inputs(vec![ios!(OSSM1Lcl6F)])
    .outputs(ios!(OSSM1Lcl, MCM2Lcl6D))
    .outputs(ios!(
        OSSAzEncoderAngle,
        OSSElEncoderAngle,
        OSSRotEncoderAngle
    ))
    .build::<Exponential>()
    .unwrap();

    let mut mount_drives_forces = Some(ios!(
        OSSAzDriveTorque(vec![0f64; 12]),
        OSSElDriveTorque(vec![0f64; 4]),
        OSSRotDriveTorque(vec![0f64; 4])
    ));

    let n_step = 30_000;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m2_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    for _ in 0..n_step {
        let mut fem_forces = vec![ios!(OSSM1Lcl6F(vec![1f64; 42]))];
        if let Some(x) = mount_drives_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        let fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap().unwrap();

        let mount_encoders = vec![
            fem_outputs[ios!(OSSElEncoderAngle)].clone(),
            fem_outputs[ios!(OSSAzEncoderAngle)].clone(),
            fem_outputs[ios!(OSSRotEncoderAngle)].clone(),
        ];
        mount_drives_forces = mnt_ctrl
            .in_step_out(Some(mount_encoders.clone()))
            .unwrap()
            .map(|mut x| {
                x.extend_from_slice(&mount_encoders);
                mnt_drives.in_step_out(Some(x))
            })
            .unwrap()
            .unwrap();

        Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)]).map(|x| m1_logs.push(x));
        Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).map(|x| m2_logs.push(x));
    }

    let mut file = File::create("tests/tests_mount_constants.pkl").unwrap();
    pkl::to_writer(&mut file, &(m1_logs, m2_logs), true).unwrap();
}
