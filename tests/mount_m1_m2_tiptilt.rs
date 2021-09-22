use dosio::{ios, Dos, IOVec, IO};
use fem::{dos::DiscreteStateSpace, FEM};
use m1_ctrl as m1;
use mount_ctrl as mount;
use nalgebra as na;
use ns_opm_im::optical_sensitivities::{from_opticals, OpticalSensitivities};
use serde_pickle as pkl;
use std::fs::File;

#[test]
fn mount_m1_m2_tiptilt_constant() {
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
        let fem = FEM::from_env().unwrap();
        println!("FEM:\n{}", fem);
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

    let opticals = from_opticals(&{ OpticalSensitivities::load().unwrap() }[3..4]);
    let mut pzt_cmd = Some(vec![ios!(TTcmd(vec![0f64; 21]))]);

    let n_fem_forces = fem_forces.len();
    let mut fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap().unwrap();
    let n_step = 10_000;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m2_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    for k in 0..n_step {
        let mut fem_forces = Vec::<IO<Vec<f64>>>::with_capacity(n_fem_forces);
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
        // Optical System
        //  - tip-tilt
        if k >= 5000 && k % 200 == 0 {
            let rbms = na::DVector::from_iterator(
                84,
                Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)])
                    .into_iter()
                    .chain(Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).into_iter())
                    .flatten(),
            );
            let stt = &opticals * rbms;
            pzt_cmd = os_tiptilt
                .in_step_out(Some(ios!(TTSP(vec![1e-6]), TTFB(stt.as_slice().to_vec()))))
                .unwrap();
        };
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
        // LOGS
        Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)]).map(|x| m1_logs.push(x));
        Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).map(|x| m2_logs.push(x));
        //        lom_logs.push(stt.as_slice().to_vec());
    }

    let segment_tiptilt: Vec<_> = m1_logs
        .iter()
        .zip(m2_logs.iter())
        .map(|(m1_rbm, m2_rbm)| {
            let rbm = na::DVector::from_iterator(
                84,
                m1_rbm.iter().cloned().chain(m2_rbm.iter().cloned()),
            );
            let stt = &opticals * rbm;
            stt.as_slice().to_vec()
        })
        .collect();
    let mut file = File::create("tests/tests_mount-m1-m2-tiptilt_constant.pkl").unwrap();
    pkl::to_writer(&mut file, &(m1_logs, m2_logs, segment_tiptilt), true).unwrap();
}
