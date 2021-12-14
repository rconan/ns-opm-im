use dosio::{io::jar, ios, Dos, IO};
use fem::{
    dos::{DiscreteStateSpace, Exponential},
    FEM,
};
use simple_logger::SimpleLogger;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    SimpleLogger::new().env().init()?;

    for k in 0..3 {
        let sampling_rate = 1e3;
        let mut fem = {
            let fem = FEM::from_env()?;
            DiscreteStateSpace::from(fem)
        }
        .sampling(sampling_rate)
        .proportional_damping(2. / 100.)
        .inputs(vec![ios!(MCM2PZTF)])
        .outputs(vec![ios!(MCM2PZTD)])
        .build::<Exponential>()?;

        let pzt_force = {
            let mut data = vec![0f64; 42];
            data[2 * k] = -1.;
            data[2 * k + 1] = 1.;
            vec![ios!(MCM2PZTF(data))]
        };

        let mut pzt_disp = None;
        for _ in 0..1_000 {
            pzt_disp =
                Option::<Vec<f64>>::from(&fem.in_step_out(Some(pzt_force.clone()))?.unwrap()[0]);
        }
        let pzt_gain: Vec<_> = pzt_disp
            .unwrap()
            .chunks(2)
            .map(|x| 1e8 * (x[1] - x[0]))
            .collect();
        println!("PZT gain: {:#?}", &pzt_gain[..3]);
    }

    Ok(())
}
