pub trait Control {
    fn step(&mut self, value: &[f64]) -> Option<Vec<f64>>;
}

#[derive(Default)]
pub struct Average {
    n_sample: usize,
    average: Vec<f64>,
    counter: usize,
}
impl Average {
    pub fn new(n_sample: usize, n_data: usize) -> Self {
        Self {
            n_sample,
            average: vec![0f64; n_data],
            counter: 0,
        }
    }
}
impl Iterator for Average {
    type Item = Vec<f64>;
    fn next(&mut self) -> Option<Self::Item> {
        self.counter += 1;
        if self.counter == self.n_sample {
            let mean: Vec<_> = self
                .average
                .iter()
                .map(|x| x / self.n_sample as f64)
                .collect();
            self.average = vec![0f64; self.average.len()];
            self.counter = 0;
            Some(mean)
        } else {
            None
        }
    }
}
impl Control for Average {
    fn step(&mut self, value: &[f64]) -> Option<Vec<f64>> {
        self.average.iter_mut().zip(value).for_each(|(a, v)| {
            *a += *v;
        });
        self.next()
    }
}

#[derive(Default)]
pub struct Integrate {
    gain: f64,
    mem: Vec<f64>,
}
impl Integrate {
    pub fn new(gain: f64, n_data: usize) -> Self {
        Self {
            gain,
            mem: vec![0f64; n_data],
        }
    }
    pub fn last(&self) -> Option<Vec<f64>> {
        Some(self.mem.clone())
    }
}
impl Control for Integrate {
    fn step(&mut self, value: &[f64]) -> Option<Vec<f64>> {
        let gain = self.gain;
        self.mem.iter_mut().zip(value).for_each(|(a, v)| {
            *a += *v * gain;
        });
        self.last()
    }
}
#[derive(Default)]
pub struct Proportional {
    gain: f64,
    mem: Vec<f64>,
}
impl Proportional {
    pub fn new(gain: f64, n_data: usize) -> Self {
        Self {
            gain,
            mem: vec![0f64; n_data],
        }
    }
    pub fn last(&self) -> Option<Vec<f64>> {
        Some(self.mem.clone())
    }
}
impl Control for Proportional {
    fn step(&mut self, value: &[f64]) -> Option<Vec<f64>> {
        let gain = self.gain;
        self.mem.iter_mut().zip(value).for_each(|(a, v)| {
            *a = *v * gain;
        });
        self.last()
    }
}

#[derive(Default)]
pub struct Delay {
    mem: Vec<Vec<f64>>,
}
impl Delay {
    pub fn new(n_sample: usize, n_data: usize) -> Self {
        Self {
            mem: vec![vec![0f64; n_data]; n_sample + 1],
        }
    }
}
impl Control for Delay {
    fn step(&mut self, value: &[f64]) -> Option<Vec<f64>> {
        self.mem[0].copy_from_slice(value);
        self.mem.rotate_right(1);
        Some(self.mem[0].clone())
    }
}

#[derive(Default)]
pub struct StateSpace2x2 {
    a: [f64; 4],
    b: [f64; 2],
    c: [f64; 2],
    d: f64,
    x: (Vec<f64>, Vec<f64>),
}
impl StateSpace2x2 {
    pub fn new(a: [f64; 4], b: [f64; 2], c: [f64; 2], d: f64, n_data: usize) -> Self {
        Self {
            a,
            b,
            c,
            d,
            x: (vec![0f64; n_data], vec![0f64; n_data]),
        }
    }
}
impl Control for StateSpace2x2 {
    fn step(&mut self, u: &[f64]) -> Option<Vec<f64>> {
        let x_iter = (&self.x.0).into_iter().zip((&self.x.1).into_iter());
        let y: Vec<_> = u
            .iter()
            .zip(x_iter.clone())
            .map(|(u, x)| self.c[0] * x.0 + self.c[1] * x.1 + self.d * u)
            .collect();
        self.x = u
            .iter()
            .zip(x_iter)
            .map(|(u, x)| {
                (
                    self.a[0] * x.0 + self.a[1] * x.1 + self.b[0] * u,
                    self.a[2] * x.0 + self.a[3] * x.1 + self.b[1] * u,
                )
            })
            .unzip();
        Some(y)
    }
}
