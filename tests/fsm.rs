use dosio::{ios, Dos, IOVec, IO};
use fem::{dos::DiscreteStateSpace, FEM};
use m1_ctrl as m1;
use mat73;
use mount_ctrl as mount;
use serde_pickle as pkl;
use skyangle::Conversion;
use std::fs::File;

#[test]
fn fsm_constant_pos() {
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
    .inputs_from(&fsm_positionner)
    .outputs(ios!(OSSM1Lcl, MCM2Lcl6D))
    .outputs(vec![ios!(MCM2SmHexD)])
    .build()
    .unwrap();

    let mut fsm_positionner_forces = Some(vec![ios!(MCM2SmHexF(vec![0f64; 84]))]);

    let n_step = 1_000;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m2_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    for _ in 0..n_step {
        // FEM
        let mut fem_forces = vec![];
        fsm_positionner_forces.as_mut().map(|x| {
            fem_forces.append(x);
        });
        let mut fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap().unwrap();
        // M2
        let data = fem_outputs.pop_these(vec![ios!(MCM2SmHexD)]).map(|mut x| {
            x.append(&mut vec![ios!(M2poscmd(vec![1e-6; 42]))]);
            x
        });
        fsm_positionner_forces = fsm_positionner.in_step_out(data).unwrap();
        // LOGS
        Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)]).map(|x| m1_logs.push(x));
        Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).map(|x| m2_logs.push(x));
    }

    let avrg_steps = 200;
    let cvgce_err = 1_f64
        - 1e6
            * m2_logs
                .iter()
                .skip(n_step - avrg_steps)
                .fold(vec![0f64; 42], |mut a, s| {
                    a.iter_mut().zip(s.iter()).for_each(|(a, s)| *a += *s);
                    a
                })
                .into_iter()
                .map(|x| x / avrg_steps as f64)
                .sum::<f64>()
            / 42_f64;
    println!("Convergence error: {}", cvgce_err);

    let mut file = File::create("tests/tests_fsm_constant-pos.pkl").unwrap();
    pkl::to_writer(&mut file, &(m1_logs, m2_logs), true).unwrap();
}

#[test]
fn fsm_constant_pzt() {
    // M2
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
    .inputs_from(&fsm_piezostack)
    .outputs(ios!(OSSM1Lcl, MCM2Lcl6D))
    .outputs(vec![ios!(MCM2PZTD)])
    .build()
    .unwrap();

    let mut fsm_piezostack_forces = Some(vec![ios!(MCM2PZTF(vec![0f64; 42]))]);

    let n_step = 1_000;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m2_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    for _ in 0..n_step {
        // FEM
        let mut fem_forces = vec![];
        fsm_piezostack_forces.as_mut().map(|x| {
            fem_forces.append(x);
        });
        let mut fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap().unwrap();
        // M2
        let data = fem_outputs.pop_these(vec![ios!(MCM2PZTD)]).map(|mut x| {
            x.append(&mut vec![ios!(TTcmd(pzt_cmd.clone()))]);
            x
        });
        fsm_piezostack_forces = fsm_piezostack.in_step_out(data).unwrap();
        // LOGS
        Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)]).map(|x| m1_logs.push(x));
        Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).map(|x| m2_logs.push(x));
    }

    let avrg_steps = 200;
    let cvgce_err = 1_f64
        - 1e6
            * m2_logs
                .iter()
                .skip(n_step - avrg_steps)
                .fold(vec![0f64; 42], |mut a, s| {
                    a.iter_mut().zip(s.iter()).for_each(|(a, s)| *a += *s);
                    a
                })
                .into_iter()
                .map(|x| x / avrg_steps as f64)
                .sum::<f64>()
            / 42_f64;
    println!("Convergence error: {}", cvgce_err);

    let mut file = File::create("tests/tests_fsm_constant-pzt.pkl").unwrap();
    pkl::to_writer(&mut file, &(m1_logs, m2_logs), true).unwrap();
}

#[test]
fn fsm_constant_tiptiltpzt() {
    // M2
    let mut fsm_piezostack = fsm::piezostack::Controller::new();
    let mut fsm_tiptilt = fsm::tiptilt::Controller::new();
    //FEM
    let sampling_rate = 1e3;
    let mut fem = {
        let fem = FEM::from_env().unwrap();
        DiscreteStateSpace::from(fem)
    }
    .sampling(sampling_rate)
    .proportional_damping(2. / 100.)
    .inputs_from(&fsm_piezostack)
    .outputs(ios!(OSSM1Lcl, MCM2Lcl6D))
    .outputs(vec![ios!(MCM2PZTD)])
    .build()
    .unwrap();

    let mut fsm_piezostack_forces = Some(vec![ios!(MCM2PZTF(vec![0f64; 42]))]);

    let n_step = 1_000;
    let mut m1_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    let mut m2_logs = Vec::<Vec<f64>>::with_capacity(n_step * 42);
    for _ in 0..n_step {
        // FEM
        let mut fem_forces = vec![];
        fsm_piezostack_forces.as_mut().map(|x| {
            fem_forces.append(x);
        });
        let mut fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap().unwrap();
        // M2
        let data = fsm_tiptilt
            .in_step_out(Some(ios!(TTSP(vec![1e-6]), TTFB(vec![0f64; 14]))))
            .unwrap()
            .map(|mut tt_cmd| {
                fem_outputs
                    .pop_these(vec![ios!(MCM2PZTD)])
                    .map(|mut pzt_fb| {
                        tt_cmd.append(&mut pzt_fb);
                        tt_cmd
                    })
            })
            .flatten();
        fsm_piezostack_forces = fsm_piezostack.in_step_out(data).unwrap();
        // LOGS
        Option::<Vec<f64>>::from(&fem_outputs[ios!(OSSM1Lcl)]).map(|x| m1_logs.push(x));
        Option::<Vec<f64>>::from(&fem_outputs[ios!(MCM2Lcl6D)]).map(|x| m2_logs.push(x));
    }

    let avrg_steps = 200;
    let cvgce_err = 1_f64
        - 1e6
            * m2_logs
                .iter()
                .skip(n_step - avrg_steps)
                .fold(vec![0f64; 42], |mut a, s| {
                    a.iter_mut().zip(s.iter()).for_each(|(a, s)| *a += *s);
                    a
                })
                .into_iter()
                .map(|x| x / avrg_steps as f64)
                .sum::<f64>()
            / 42_f64;
    println!("Convergence error: {}", cvgce_err);

    let mut file = File::create("tests/tests_fsm_constant-tiptiltpzt.pkl").unwrap();
    pkl::to_writer(&mut file, &(m1_logs, m2_logs), true).unwrap();
}
