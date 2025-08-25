pub mod models {
    use num_complex::Complex64;

    pub trait CharacteristicFunction {
        fn phi(&self, u: Complex64, t: f64, s0: f64, r: f64, q: f64) -> Complex64;
    }

    #[derive(Clone, Copy, Debug)]
    pub struct BlackScholes {
        pub sigma: f64,
    }

    impl CharacteristicFunction for BlackScholes {
        fn phi(&self, u: Complex64, t: f64, s0: f64, r: f64, q: f64) -> Complex64 {
            let sigma2 = self.sigma * self.sigma;
            let iu = Complex64::new(0.0, 1.0) * u;
            let drift = (r - q - 0.5 * sigma2) * t + s0.ln();
            (iu * drift - 0.5 * sigma2 * u * u * t).exp()
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct Heston {
        pub kappa: f64,
        pub theta: f64,
        pub vol_of_vol: f64,
        pub rho: f64,
        pub v0: f64,
    }

    impl CharacteristicFunction for Heston {
        fn phi(&self, u: Complex64, t: f64, s0: f64, r: f64, q: f64) -> Complex64 {
            let i = Complex64::new(0.0, 1.0);
            let iu = i * u;
            let sigma = self.vol_of_vol;
            let kappa = self.kappa;
            let theta = self.theta;
            let rho = self.rho;
            let v0 = self.v0;

            let d = ((rho * sigma * iu - kappa).powi(2) + sigma * sigma * (iu + u * u)).sqrt();
            let g = (kappa - rho * sigma * iu - d) / (kappa - rho * sigma * iu + d);

            let exp_neg_d_t = (-d * t).exp();
            let one = Complex64::new(1.0, 0.0);

            let c = (kappa * theta / (sigma * sigma))
                * ((kappa - rho * sigma * iu - d) * t - 2.0 * ((one - g * exp_neg_d_t) / (one - g)).ln());
            let d_term = ((kappa - rho * sigma * iu - d) / (sigma * sigma)) * ((one - exp_neg_d_t) / (one - g * exp_neg_d_t));

            (iu * (s0.ln() + (r - q) * t) + v0 * d_term + c).exp()
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct VarianceGamma {
        pub sigma: f64,
        pub nu: f64,
        pub theta: f64,
    }

    impl CharacteristicFunction for VarianceGamma {
        fn phi(&self, u: Complex64, t: f64, s0: f64, r: f64, q: f64) -> Complex64 {
            let i = Complex64::new(0.0, 1.0);
            let iu = i * u;
            let sigma2 = self.sigma * self.sigma;
            let nu = self.nu;
            let theta = self.theta;

            let omega = -((1.0 - theta * nu - 0.5 * sigma2 * nu).ln()) / nu;
            let drift = s0.ln() + (r - q + omega) * t;

            let base = Complex64::new(1.0, 0.0) - i * theta * nu * u + Complex64::new(0.5 * sigma2 * nu, 0.0) * u * u;
            (iu * drift).exp() * base.powf(-t / nu)
        }
    }
}

pub mod fft_pricer {
    use num_complex::Complex64;
    use rustfft::FftPlanner;
    use std::f64::consts::PI;

    use crate::models::CharacteristicFunction;

    pub struct CarrMadanParams {
        pub alpha: f64,
        pub eta: f64,
        pub n: usize, // grid points is N = 2^n
        pub beta: f64, // log-strike shift, commonly ln(K)
    }

    pub struct GridResult {
        pub k: Vec<f64>,           // log-strikes
        pub call_prices: Vec<f64>, // call prices on the grid
    }

    pub fn price_calls_grid<M: CharacteristicFunction>(
        model: &M,
        s0: f64,
        r: f64,
        q: f64,
        t: f64,
        params: CarrMadanParams,
    ) -> GridResult {
        let n_points = 1usize << params.n;
        let eta = params.eta;
        let alpha = params.alpha;
        let beta = params.beta;

        let lambda = 2.0 * PI / ((n_points as f64) * eta);

        let mut y: Vec<Complex64> = Vec::with_capacity(n_points);
        let i = Complex64::new(0.0, 1.0);
        let discount = (-r * t).exp();

        for j in 0..n_points {
            let u_j = (j as f64) * eta;
            let u_shifted = Complex64::new(u_j, -(alpha + 1.0)); // u - i(alpha+1)

            let phi_val = model.phi(u_shifted, t, s0, r, q);
            let denom = (alpha * alpha + alpha - u_j * u_j) + i * ((2.0 * alpha + 1.0) * u_j);
            let psi = discount * phi_val / denom;

            let weight = if j == 0 { 1.0 } else if j % 2 == 1 { 4.0 } else { 2.0 };
            let weight = weight / 3.0;

            let factor = (i * (beta * u_j)).exp();
            let val = psi * factor * (eta * weight);
            y.push(val);
        }

        let mut planner = FftPlanner::<f64>::new();
        let fft = planner.plan_fft_forward(n_points);
        let mut y_fft = y;
        fft.process(&mut y_fft);

        let mut k_grid: Vec<f64> = Vec::with_capacity(n_points);
        let mut call_prices: Vec<f64> = Vec::with_capacity(n_points);
        for m in 0..n_points {
            let k_m = -beta + (m as f64) * lambda;
            let c_m = (-(alpha * k_m)).exp() * (y_fft[m].re) / PI;
            k_grid.push(k_m);
            call_prices.push(c_m);
        }

        GridResult { k: k_grid, call_prices }
    }

    pub fn price_put_at_strike<M: CharacteristicFunction>(
        model: &M,
        s0: f64,
        r: f64,
        q: f64,
        t: f64,
        params: CarrMadanParams,
        k_strike: f64,
    ) -> f64 {
        let grid = price_calls_grid(model, s0, r, q, t, params);
        let target_k = k_strike.ln();
        let mut best_idx = 0usize;
        let mut best_err = f64::INFINITY;
        for (idx, &k) in grid.k.iter().enumerate() {
            let err = (k - target_k).abs();
            if err < best_err {
                best_err = err;
                best_idx = idx;
            }
        }
        let call_price = grid.call_prices[best_idx];
        // Put from parity: P = C - S0 e^{-qT} + K e^{-rT}
        let put_price = call_price - s0 * (-q * t).exp() + k_strike * (-r * t).exp();
        put_price
    }
} 