use dosio::{ios, Dos, IOVec, IO};
use fem::{
    dos::{DiscreteModalSolver, DiscreteStateSpace, Exponential},
    FEM,
};
use geotrans::{Quaternion, Vector};
use m1_ctrl as m1;
use mount_ctrl as mount;
use nalgebra as na;
use serde_pickle as pkl;
use skyangle::Conversion;
use std::fs::File;

#[test]
fn mount_m1_zeros() {
    // MOUNT
    let mut mnt_drives = mount::drives::Controller::new();
    let mut mnt_ctrl = mount::controller::Controller::new();
    // M1
    let mut m1_hardpoints = m1::hp_dynamics::Controller::new();
    let mut m1_load_cells = m1::hp_load_cells::Controller::new();
    let mut m1_actuators = m1::actuators::M1ForceLoops::new();
    //FEM
    let sampling_rate = 1e3;
    let mut fem = {
        let fem = FEM::from_env().unwrap();
        DiscreteStateSpace::from(fem)
            .sampling(sampling_rate)
            .proportional_damping(2. / 100.)
            .inputs_from(&mnt_drives)
            .inputs_from(&m1_hardpoints)
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
            .build()
            .unwrap()
    };

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

    for k in 0..3 {
        let mut fem_forces = vec![];
        if let Some(x) = mount_drives_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        if let Some(x) = m1_hardpoints_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        if let Some(x) = m1_actuators_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        let mut fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap().unwrap();

        let mount_encoders = <Vec<IO<Vec<f64>>> as IOVec>::pop_these(
            &mut fem_outputs,
            ios!(OSSElEncoderAngle, OSSAzEncoderAngle, OSSRotEncoderAngle),
        );
        mount_drives_forces = mnt_ctrl
            .in_step_out(mount_encoders.clone())
            .unwrap()
            .map(|mut x| {
                x.extend_from_slice(mount_encoders.as_ref().unwrap());
                mnt_drives.in_step_out(Some(x))
            })
            .unwrap()
            .unwrap();

        let data =
            <Vec<IO<Vec<f64>>> as IOVec>::pop_these(&mut fem_outputs, vec![ios!(OSSHardpointD)])
                .map(|mut x| {
                    x.extend_from_slice(m1_hardpoints_forces.as_ref().unwrap());
                    x
                });
        let mut m1_hp_lc = m1_load_cells.in_step_out(data).unwrap();
        if k % 10 == 0 {
            if let Some(ref mut m1_hp_lc) = m1_hp_lc {
                m1_hp_lc.append(&mut m1_bending_modes.clone());
            }
            m1_actuators_forces = m1_actuators.in_step_out(m1_hp_lc).unwrap();
        };
        m1_hardpoints_forces = m1_hardpoints
            .in_step_out(Some(vec![ios!(M1RBMcmd(vec![0f64; 42]))]))
            .unwrap();
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
fn mount_m1_constants() {
    // MOUNT
    let mut mnt_drives = mount::drives::Controller::new();
    let mut mnt_ctrl = mount::controller::Controller::new();
    // M1
    let mut m1_hardpoints = m1::hp_dynamics::Controller::new();
    let mut m1_load_cells = m1::hp_load_cells::Controller::new();
    let mut m1_actuators = m1::actuators::M1ForceLoops::new();
    //FEM
    let sampling_rate = 1e3;
    let mut fem = {
        let fem = FEM::from_env().unwrap();
        DiscreteStateSpace::from(fem)
    }
    .sampling(sampling_rate)
    .proportional_damping(2. / 100.)
    .inputs_from(&mnt_drives)
    .inputs_from(&m1_hardpoints)
    .inputs(vec![ios!(OSSM1Lcl6F)])
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
    .build()
    .unwrap();

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

    let n_step = 30_000;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m2_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    for k in 0..n_step {
        // FEM
        let mut fem_forces = vec![ios!(OSSM1Lcl6F(vec![1f64; 42]))];
        if let Some(x) = mount_drives_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        if let Some(x) = m1_hardpoints_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        if let Some(x) = m1_actuators_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        let mut fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap().unwrap();
        // MOUNT
        let mount_encoders = <Vec<IO<Vec<f64>>> as IOVec>::pop_these(
            &mut fem_outputs,
            ios!(OSSElEncoderAngle, OSSAzEncoderAngle, OSSRotEncoderAngle),
        );
        mount_drives_forces = mnt_ctrl
            .in_step_out(mount_encoders.clone())
            .unwrap()
            .map(|mut x| {
                x.extend_from_slice(mount_encoders.as_ref().unwrap());
                mnt_drives.in_step_out(Some(x))
            })
            .unwrap()
            .unwrap();
        // M1
        let data =
            <Vec<IO<Vec<f64>>> as IOVec>::pop_these(&mut fem_outputs, vec![ios!(OSSHardpointD)])
                .map(|mut x| {
                    x.extend_from_slice(m1_hardpoints_forces.as_ref().unwrap());
                    x
                });
        let mut m1_hp_lc = m1_load_cells.in_step_out(data).unwrap();
        if k % 10 == 0 {
            if let Some(ref mut m1_hp_lc) = m1_hp_lc {
                m1_hp_lc.append(&mut m1_bending_modes.clone());
            }
            m1_actuators_forces = m1_actuators.in_step_out(m1_hp_lc).unwrap();
        };
        m1_hardpoints_forces = m1_hardpoints
            .in_step_out(Some(vec![ios!(M1RBMcmd(vec![0f64; 42]))]))
            .unwrap();
        // LOGS
        Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)]).map(|x| m1_logs.push(x));
        Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).map(|x| m2_logs.push(x));
    }

    let mut file = File::create("tests/tests_mount-m1_constants.pkl").unwrap();
    pkl::to_writer(&mut file, &(m1_logs, m2_logs), true).unwrap();
}

#[test]
fn mount_m1_rbms() {
    // MOUNT
    let mut mnt_drives = mount::drives::Controller::new();
    let mut mnt_ctrl = mount::controller::Controller::new();
    // M1
    let mut m1_hardpoints = m1::hp_dynamics::Controller::new();
    let mut m1_load_cells = m1::hp_load_cells::Controller::new();
    let mut m1_actuators = m1::actuators::M1ForceLoops::new();
    //FEM
    let sampling_rate = 1e3;
    let (m1_segments_surface_nodes, mut fem): (Vec<Vec<f64>>, DiscreteModalSolver<Exponential>) = {
        let fem = FEM::from_env().unwrap();
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
                .inputs_from(&mnt_drives)
                .inputs_from(&m1_hardpoints)
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
                .outputs(ios!(
                    M1Segment1AxialD,
                    M1Segment2AxialD,
                    M1Segment3AxialD,
                    M1Segment4AxialD,
                    M1Segment5AxialD,
                    M1Segment6AxialD,
                    M1Segment7AxialD
                ))
                .build()
                .unwrap(),
        )
    };

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
        M1S1BMcmd(vec![0f64; 27]),
        M1S2BMcmd(vec![0f64; 27]),
        M1S3BMcmd(vec![0f64; 27]),
        M1S4BMcmd(vec![0f64; 27]),
        M1S5BMcmd(vec![0f64; 27]),
        M1S6BMcmd(vec![0f64; 27]),
        M1S7BMcmd(vec![0f64; 27])
    );
    let mut m1_hardpoints_forces = Some(vec![ios!(OSSHarpointDeltaF(vec![0f64; 42]))]);

    let n_step = 5_000;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m2_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m1_shape = vec![];
    for k in 0..n_step {
        // FEM
        let mut fem_forces = vec![];
        if let Some(x) = mount_drives_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        if let Some(x) = m1_hardpoints_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        if let Some(x) = m1_actuators_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        let mut fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap().unwrap();
        // MOUNT
        let mount_encoders = fem_outputs.pop_these(ios!(
            OSSElEncoderAngle,
            OSSAzEncoderAngle,
            OSSRotEncoderAngle
        ));
        mount_drives_forces = mnt_ctrl
            .in_step_out(mount_encoders.clone())
            .unwrap()
            .map(|mut x| {
                x.extend_from_slice(mount_encoders.as_ref().unwrap());
                mnt_drives.in_step_out(Some(x))
            })
            .unwrap()
            .unwrap();
        // M1
        let data = fem_outputs
            .pop_these(vec![ios!(OSSHardpointD)])
            .map(|mut x| {
                x.extend_from_slice(m1_hardpoints_forces.as_ref().unwrap());
                x
            });
        let mut m1_hp_lc = m1_load_cells.in_step_out(data).unwrap();
        if k % 10 == 0 {
            if let Some(ref mut m1_hp_lc) = m1_hp_lc {
                m1_hp_lc.append(&mut m1_bending_modes.clone());
            }
            m1_actuators_forces = m1_actuators.in_step_out(m1_hp_lc).unwrap();
        };

        m1_hardpoints_forces = m1_hardpoints
            .in_step_out(Some(if k < 2000 {
                vec![ios!(M1RBMcmd(vec![0f64; 42]))]
            } else {
                let mut rbms = vec![0f64; 42];
                rbms[3] = 1e-6;
                vec![ios!(M1RBMcmd(rbms))]
            }))
            .unwrap();
        // LOGS
        Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)]).map(|x| m1_logs.push(x));
        Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).map(|x| m2_logs.push(x));
        m1_shape = ios!(
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
        .zip(m1_segments_surface_nodes.iter())
        .zip(m1_logs.last().unwrap().chunks(6))
        .map(|((surface, nodes), rbm)| {
            let (t_xyz, r_xyz) = rbm.split_at(3);
            let rxyz_surface = {
                // 3D rotation of mirror surface nodes
                let q = Quaternion::unit(r_xyz[2], Vector::k())
                    * Quaternion::unit(r_xyz[1], Vector::j())
                    * Quaternion::unit(r_xyz[0], Vector::i());
                let trans_nodes: Vec<f64> = nodes
                    .chunks(3)
                    .flat_map(|x| {
                        let p: Quaternion = Vector::from(x).into();
                        let w: Quaternion = &q * p * &q.complex_conjugate();
                        let vv: Vec<f64> = (*Vector::from(w.vector_as_slice())).into();
                        vv
                    })
                    .collect();
                trans_nodes
                    .chunks(3)
                    .map(|x| x[2])
                    .zip(nodes.chunks(3).map(|x| x[2]))
                    .map(|(z_rbm, z)| z_rbm - z)
                    .collect::<Vec<f64>>()
            };
            // Removing Rx, Ry, Rz and Tz
            surface
                .iter()
                .zip(rxyz_surface.iter())
                .map(|(a, r)| a - r - t_xyz[2])
                .collect::<Vec<f64>>()
        })
        .collect();
    }

    let mut file = File::create("tests/tests_mount-m1_rbms.pkl").unwrap();
    pkl::to_writer(&mut file, &(m1_logs, m2_logs, m1_shape), true).unwrap();
}
