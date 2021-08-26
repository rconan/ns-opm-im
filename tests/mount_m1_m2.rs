use dosio::{ios, Dos, DosVec, IOVec, IO};
use fem::{dos::DiscreteStateSpace, FEM};
use m1_ctrl as m1;
use mount_ctrl as mount;
use serde_pickle as pkl;
use skyangle::Conversion;
use std::fs::File;

#[test]
fn mount_m1_m2_zeros() {
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
    .inputs(ios!(
        M1ActuatorsSegment1,
        M1ActuatorsSegment2,
        M1ActuatorsSegment3,
        M1ActuatorsSegment4,
        M1ActuatorsSegment5,
        M1ActuatorsSegment6,
        M1ActuatorsSegment7
    ))
    .inputs_from(&fsm_positionner)
    .inputs_from(&fsm_piezostack)
    .outputs(ios!(OSSM1Lcl, MCM2Lcl6D))
    .outputs(ios!(
        OSSAzEncoderAngle,
        OSSElEncoderAngle,
        OSSRotEncoderAngle
    ))
    .outputs(vec![ios!(OSSHardpointD)])
    .outputs(ios!(MCM2SmHexD, MCM2PZTD))
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
        M1S1BMcmd(vec![0f64; 27]),
        M1S2BMcmd(vec![0f64; 27]),
        M1S3BMcmd(vec![0f64; 27]),
        M1S4BMcmd(vec![0f64; 27]),
        M1S5BMcmd(vec![0f64; 27]),
        M1S6BMcmd(vec![0f64; 27]),
        M1S7BMcmd(vec![0f64; 27])
    );
    let mut m1_hardpoints_forces = Some(vec![ios!(OSSHarpointDeltaF(vec![0f64; 42]))]);
    let mut fsm_positionner_forces = Some(vec![ios!(MCM2SmHexF(vec![0f64; 84]))]);
    let mut fsm_piezostack_forces = Some(vec![ios!(MCM2PZTF(vec![0f64; 42]))]);

    for k in 0..3 {
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
        if let Some(x) = fsm_positionner_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        if let Some(x) = fsm_piezostack_forces.as_mut() {
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
        // M2
        let data =
            <Vec<IO<Vec<f64>>> as IOVec>::pop_these(&mut fem_outputs, vec![ios!(MCM2SmHexD)]).map(
                |mut x| {
                    x.append(&mut vec![ios!(M2poscmd(vec![0f64; 42]))]);
                    x
                },
            );
        fsm_positionner_forces = fsm_positionner.in_step_out(data).unwrap();
        let data = <Vec<IO<Vec<f64>>> as IOVec>::pop_these(&mut fem_outputs, vec![ios!(MCM2PZTD)])
            .map(|mut x| {
                x.append(&mut vec![ios!(TTcmd(vec![0f64; 21]))]);
                x
            });
        fsm_piezostack_forces = fsm_piezostack.in_step_out(data).unwrap();
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
fn mount_m1_m2_constant_windloads() {
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
    .inputs_from(&fsm_positionner)
    .inputs_from(&fsm_piezostack)
    .outputs(ios!(OSSM1Lcl, MCM2Lcl6D))
    .outputs(ios!(
        OSSAzEncoderAngle,
        OSSElEncoderAngle,
        OSSRotEncoderAngle
    ))
    .outputs(vec![ios!(OSSHardpointD)])
    .outputs(ios!(MCM2SmHexD, MCM2PZTD))
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
        M1S1BMcmd(vec![0f64; 27]),
        M1S2BMcmd(vec![0f64; 27]),
        M1S3BMcmd(vec![0f64; 27]),
        M1S4BMcmd(vec![0f64; 27]),
        M1S5BMcmd(vec![0f64; 27]),
        M1S6BMcmd(vec![0f64; 27]),
        M1S7BMcmd(vec![0f64; 27])
    );
    let mut m1_hardpoints_forces = Some(vec![ios!(OSSHarpointDeltaF(vec![0f64; 42]))]);
    let mut fsm_positionner_forces = Some(vec![ios!(MCM2SmHexF(vec![0f64; 84]))]);
    let mut fsm_piezostack_forces = Some(vec![ios!(MCM2PZTF(vec![0f64; 42]))]);

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
        if let Some(x) = fsm_positionner_forces.as_mut() {
            fem_forces.append(&mut x.clone());
        }
        if let Some(x) = fsm_piezostack_forces.as_mut() {
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
        // M2
        let data =
            <Vec<IO<Vec<f64>>> as IOVec>::pop_these(&mut fem_outputs, vec![ios!(MCM2SmHexD)]).map(
                |mut x| {
                    x.append(&mut vec![ios!(M2poscmd(vec![0f64; 42]))]);
                    x
                },
            );
        fsm_positionner_forces = fsm_positionner.in_step_out(data).unwrap();
        let data = <Vec<IO<Vec<f64>>> as IOVec>::pop_these(&mut fem_outputs, vec![ios!(MCM2PZTD)])
            .map(|mut x| {
                x.append(&mut vec![ios!(TTcmd(vec![0f64; 21]))]);
                x
            });
        fsm_piezostack_forces = fsm_piezostack.in_step_out(data).unwrap();
        // LOGS
        Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)]).map(|x| m1_logs.push(x));
        Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).map(|x| m2_logs.push(x));
    }

    let mut file = File::create("tests/tests_mount-m1-m2_constants.pkl").unwrap();
    pkl::to_writer(&mut file, &(m1_logs, m2_logs), true).unwrap();
}

#[test]
fn mount_m1_m2_constant_fsmpos() {
    // MOUNT
    let mut mnt_drives = mount::drives::Controller::new();
    let mut mnt_ctrl = mount::controller::Controller::new();
    // M1
    let mut m1_hardpoints = m1::hp_dynamics::Controller::new();
    let mut m1_load_cells = m1::hp_load_cells::Controller::new();
    let mut m1_actuators = m1::actuators::M1ForceLoops::new();
    // M2
    let mut fsm_positionner = fsm::positionner::Controller::new();
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
    .inputs(ios!(
        M1ActuatorsSegment1,
        M1ActuatorsSegment2,
        M1ActuatorsSegment3,
        M1ActuatorsSegment4,
        M1ActuatorsSegment5,
        M1ActuatorsSegment6,
        M1ActuatorsSegment7
    ))
    .inputs_from(&fsm_positionner)
    .outputs(ios!(OSSM1Lcl, MCM2Lcl6D))
    .outputs(ios!(
        OSSAzEncoderAngle,
        OSSElEncoderAngle,
        OSSRotEncoderAngle
    ))
    .outputs(vec![ios!(OSSHardpointD)])
    .outputs(vec![ios!(MCM2SmHexD)])
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
        M1S1BMcmd(vec![0f64; 27]),
        M1S2BMcmd(vec![0f64; 27]),
        M1S3BMcmd(vec![0f64; 27]),
        M1S4BMcmd(vec![0f64; 27]),
        M1S5BMcmd(vec![0f64; 27]),
        M1S6BMcmd(vec![0f64; 27]),
        M1S7BMcmd(vec![0f64; 27])
    );
    let mut m1_hardpoints_forces = Some(vec![ios!(OSSHarpointDeltaF(vec![0f64; 42]))]);
    let mut fsm_positionner_forces = Some(vec![ios!(MCM2SmHexF(vec![0f64; 84]))]);

    let n_step = 30_000;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m2_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
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
        if let Some(x) = fsm_positionner_forces.as_mut() {
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
        // M2
        let data =
            <Vec<IO<Vec<f64>>> as IOVec>::pop_these(&mut fem_outputs, vec![ios!(MCM2SmHexD)]).map(
                |mut x| {
                    x.append(&mut vec![ios!(M2poscmd(vec![1e-6; 42]))]);
                    x
                },
            );
        fsm_positionner_forces = fsm_positionner.in_step_out(data).unwrap();
        // LOGS
        Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)]).map(|x| m1_logs.push(x));
        Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).map(|x| m2_logs.push(x));
    }

    let mut file = File::create("tests/tests_mount-m1-m2_constant-fsmpos.pkl").unwrap();
    pkl::to_writer(&mut file, &(m1_logs, m2_logs), true).unwrap();
}

#[test]
fn mount_m1_m2_constant_fsmpospzt() {
    // MOUNT
    let mut mnt_drives = mount::drives::Controller::new();
    let mut mnt_ctrl = mount::controller::Controller::new();
    // M1
    let mut m1_hardpoints = m1::hp_dynamics::Controller::new();
    let mut m1_load_cells = m1::hp_load_cells::Controller::new();
    let mut m1_actuators = m1::actuators::M1ForceLoops::new();
    // M2
    let mut fsm_positionner = fsm::positionner::Controller::new();
    let pzt_cmd: Vec<_> = {
        let tt2pzt: Vec<f64> = {
            let mat_file = mat73::File::new("data/tt2pzt.mat").unwrap();
            mat_file.array("tt2pzt").unwrap().into()
        };
        tt2pzt.chunks(21).fold(vec![0f64; 21], |mut a, c| {
            a.iter_mut().zip(c.iter()).for_each(|(a, c)| {
                *a += *c * 1e-6;
            });
            a
        })
    };
    let mut fsm_piezostack = fsm::piezostack::Controller::new();
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
    .inputs(ios!(
        M1ActuatorsSegment1,
        M1ActuatorsSegment2,
        M1ActuatorsSegment3,
        M1ActuatorsSegment4,
        M1ActuatorsSegment5,
        M1ActuatorsSegment6,
        M1ActuatorsSegment7
    ))
    .inputs_from(&fsm_positionner)
    .inputs_from(&fsm_piezostack)
    .outputs(ios!(OSSM1Lcl, MCM2Lcl6D))
    .outputs(ios!(
        OSSAzEncoderAngle,
        OSSElEncoderAngle,
        OSSRotEncoderAngle
    ))
    .outputs(vec![ios!(OSSHardpointD)])
    .outputs(ios!(MCM2SmHexD, MCM2PZTD))
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
        M1S1BMcmd(vec![0f64; 27]),
        M1S2BMcmd(vec![0f64; 27]),
        M1S3BMcmd(vec![0f64; 27]),
        M1S4BMcmd(vec![0f64; 27]),
        M1S5BMcmd(vec![0f64; 27]),
        M1S6BMcmd(vec![0f64; 27]),
        M1S7BMcmd(vec![0f64; 27])
    );
    let mut m1_hardpoints_forces = Some(vec![ios!(OSSHarpointDeltaF(vec![0f64; 42]))]);
    let mut fsm_positionner_forces = Some(vec![ios!(MCM2SmHexF(vec![0f64; 84]))]);
    let mut fsm_piezostack_forces = Some(vec![ios!(MCM2PZTF(vec![0f64; 42]))]);

    let n_step = 30_000;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m2_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    for k in 0..n_step {
        // FEM
        let mut fem_forces = vec![];
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
        let mut fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap().unwrap();
        // MOUNT
        mount_drives_forces = DosVec::<IO<Vec<f64>>>::pop_these(
            &mut fem_outputs,
            ios!(OSSElEncoderAngle, OSSAzEncoderAngle, OSSRotEncoderAngle),
        )
        .and_then(|mut mnt_encdr| {
            mnt_ctrl
                .in_step_out(Some(mnt_encdr.clone()))
                .unwrap()
                .and_then(|mut mnt_cmd| {
                    mnt_cmd.append(&mut mnt_encdr);
                    mnt_drives.in_step_out(Some(mnt_cmd)).unwrap()
                })
        });
        // M1
        let mut m1_hp_lc =
            DosVec::<IO<Vec<f64>>>::pop_these(&mut fem_outputs, vec![ios!(OSSHardpointD)])
                .and_then(|mut hp_d| {
                    m1_hardpoints_forces.as_mut().and_then(|hp_f| {
                        hp_d.append(hp_f);
                        m1_load_cells.in_step_out(Some(hp_d)).unwrap()
                    })
                });
        if k % 10 == 0 {
            m1_actuators_forces = m1_hp_lc.and_then(|mut hp_lc| {
                hp_lc.append(&mut m1_bending_modes.clone());
                m1_actuators.in_step_out(Some(hp_lc)).unwrap()
            });
        };
        m1_hardpoints_forces = m1_hardpoints
            .in_step_out(Some(vec![ios!(M1RBMcmd(vec![0f64; 42]))]))
            .unwrap();
        // M2
        //  - positioner
        fsm_positionner_forces =
            DosVec::<IO<Vec<f64>>>::pop_these(&mut fem_outputs, vec![ios!(MCM2SmHexD)]).and_then(
                |mut hex_d| {
                    hex_d.append(&mut vec![ios!(M2poscmd(vec![0f64; 42]))]);
                    fsm_positionner.in_step_out(Some(hex_d)).unwrap()
                },
            );
        //  - piezostack
        fsm_piezostack_forces =
            DosVec::<IO<Vec<f64>>>::pop_these(&mut fem_outputs, vec![ios!(MCM2PZTD)]).and_then(
                |mut pzt_d| {
                    pzt_d.append(&mut vec![ios!(TTcmd(pzt_cmd.clone()))]);
                    fsm_piezostack.in_step_out(Some(pzt_d)).unwrap()
                },
            );
        // LOGS
        Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)]).map(|x| m1_logs.push(x));
        Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).map(|x| m2_logs.push(x));
    }

    let mut file = File::create("tests/tests_mount-m1-m2_constant-fsmpospzt.pkl").unwrap();
    pkl::to_writer(&mut file, &(m1_logs, m2_logs), true).unwrap();
}
