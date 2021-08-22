use criterion::{criterion_group, criterion_main, Criterion};
use crseo::{
    dos::GmtOpticalSensorModel,
    shackhartmann::{Diffractive, Geometric, WavefrontSensorBuilder},
    Builder, Calibration, OpticalSensitivities, ShackHartmann, SH48,
};
use dosio::{ios, Dos};
use fem::{
    dos::{DiscreteModalSolver, DiscreteStateSpace, Exponential},
    FEM,
};
use std::path::{Path, PathBuf};

fn load_fem(path: PathBuf) -> DiscreteModalSolver<Exponential> {
    let sampling_rate = 1e3;
    {
        let fem = FEM::from_pickle(path.join("modal_state_space_model_2ndOrder.73.pkl")).unwrap();
        DiscreteStateSpace::from(fem)
    }
    .sampling(sampling_rate)
    .proportional_damping(2. / 100.)
    .inputs(ios!(
        M1ActuatorsSegment1,
        M1ActuatorsSegment2,
        M1ActuatorsSegment3,
        M1ActuatorsSegment4,
        M1ActuatorsSegment5,
        M1ActuatorsSegment6,
        M1ActuatorsSegment7
    ))
    .inputs(ios!(OSSAzDriveTorque, OSSElDriveTorque, OSSRotDriveTorque))
    .inputs(vec![ios!(OSSTruss6F)])
    .outputs(ios!(
        M1Segment1AxialD,
        M1Segment2AxialD,
        M1Segment3AxialD,
        M1Segment4AxialD,
        M1Segment5AxialD,
        M1Segment6AxialD,
        M1Segment7AxialD
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
}

fn fem_benchmark(c: &mut Criterion) {
    let mut fem = load_fem(Path::new("data").join("20210802_0755_MT_mount_v202104_FSM/xyz/500Hz/"));
    c.bench_function("FEM 20210802_0755_MT_mount_v202104_FSM/xyz/500Hz/", |b| {
        b.iter(|| {
            let fem_forces = ios!(
                OSSTruss6F(vec![0f64; 6]),
                OSSAzDriveTorque(vec![0f64; 12]),
                OSSElDriveTorque(vec![0f64; 4]),
                OSSRotDriveTorque(vec![0f64; 4]),
                M1ActuatorsSegment1(vec![0f64; 335]),
                M1ActuatorsSegment2(vec![0f64; 335]),
                M1ActuatorsSegment3(vec![0f64; 335]),
                M1ActuatorsSegment4(vec![0f64; 335]),
                M1ActuatorsSegment5(vec![0f64; 335]),
                M1ActuatorsSegment6(vec![0f64; 335]),
                M1ActuatorsSegment7(vec![0f64; 306])
            );
            let _fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap();
        })
    });
}
fn fem_75Hz_benchmark(c: &mut Criterion) {
    let mut fem = load_fem(Path::new("data").join("20210802_0755_MT_mount_v202104_FSM/xyz/500Hz/"));
    c.bench_function("FEM 20210802_0755_MT_mount_v202104_FSM/xyz/75Hz/", |b| {
        b.iter(|| {
            let fem_forces = ios!(
                OSSTruss6F(vec![0f64; 6]),
                OSSAzDriveTorque(vec![0f64; 12]),
                OSSElDriveTorque(vec![0f64; 4]),
                OSSRotDriveTorque(vec![0f64; 4]),
                M1ActuatorsSegment1(vec![0f64; 335]),
                M1ActuatorsSegment2(vec![0f64; 335]),
                M1ActuatorsSegment3(vec![0f64; 335]),
                M1ActuatorsSegment4(vec![0f64; 335]),
                M1ActuatorsSegment5(vec![0f64; 335]),
                M1ActuatorsSegment6(vec![0f64; 335]),
                M1ActuatorsSegment7(vec![0f64; 306])
            );
            let _fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap();
        })
    });
}
fn fem0_benchmark(c: &mut Criterion) {
    let mut fem = load_fem(Path::new("data").join("20210802_0755_MT_mount_v202104_FSM/xyz/500Hz/"));
    c.bench_function("FEM 20210802_0755_MT_mount_v202104_FSM/", |b| {
        b.iter(|| {
            let fem_forces = ios!(
                OSSTruss6F(vec![0f64; 6]),
                OSSAzDriveTorque(vec![0f64; 12]),
                OSSElDriveTorque(vec![0f64; 4]),
                OSSRotDriveTorque(vec![0f64; 4]),
                M1ActuatorsSegment1(vec![0f64; 335]),
                M1ActuatorsSegment2(vec![0f64; 335]),
                M1ActuatorsSegment3(vec![0f64; 335]),
                M1ActuatorsSegment4(vec![0f64; 335]),
                M1ActuatorsSegment5(vec![0f64; 335]),
                M1ActuatorsSegment6(vec![0f64; 335]),
                M1ActuatorsSegment7(vec![0f64; 306])
            );
            let _fem_outputs = fem.in_step_out(Some(fem_forces)).unwrap();
        })
    });
}
fn wfsing(c: &mut Criterion) {
    type WFS_TYPE = Diffractive;
    let mut gosm = GmtOpticalSensorModel::<ShackHartmann<WFS_TYPE>, SH48<WFS_TYPE>>::new()
        .sensor(SH48::<WFS_TYPE>::new())
        .atmosphere(crseo::ATMOSPHERE::new().ray_tracing(
            26.,
            520,
            0.,
            25.,
            Some("ns-opm-im_atm.bin".to_string()),
            Some(8),
        ))
        .build()
        .unwrap();
    c.bench_function("SH48", |b| {
        b.iter(|| gosm.in_step_out(Some(vec![ios!(MCM2Lcl6D(vec![0f64; 42]))])))
    });
}

criterion_group!(
    benches,
    fem0_benchmark,
    fem_75Hz_benchmark,
    fem_benchmark,
    wfsing
);
criterion_main!(benches);
