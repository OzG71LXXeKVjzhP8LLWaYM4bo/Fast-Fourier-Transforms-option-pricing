use clap::Parser;
use fft::fft_pricer::{CarrMadanParams, price_put_at_strike};
use fft::models::VarianceGamma;

#[derive(Parser, Debug)]
#[command(name = "vg", about = "Variance-Gamma put pricing via Carr–Madan FFT")] 
struct Args {
    #[arg(long, default_value = "100.0")] s0: f64,
    #[arg(long, default_value = "80.0")] k: f64,
    #[arg(long, default_value = "0.055")] r: f64,
    #[arg(long, default_value = "0.03")] q: f64,
    #[arg(long, default_value = "1.0")] t: f64,
    #[arg(long, default_value = "0.3")] sigma: f64,
    #[arg(long, default_value = "0.5")] nu: f64,
    #[arg(long, default_value = "-0.4")] theta: f64,
}

fn main() {
    let args = Args::parse();
    let model = VarianceGamma { sigma: args.sigma, nu: args.nu, theta: args.theta };

    let alpha_values = vec![1.01, 1.25, 1.50, 1.75, 2.00, 5.00];
    let eta_values = vec![0.10, 0.25];
    let n_values = vec![6usize, 10usize];

    println!("Model = VG");
    println!("eta\tN\talpha\tput");

    for &eta in &eta_values {
        for &n in &n_values {
            for &alpha in &alpha_values {
                let params = CarrMadanParams { alpha, eta, n, beta: args.k.ln() };
                let put = price_put_at_strike(&model, args.s0, args.r, args.q, args.t, params, args.k);
                println!("{:.2}\t2^{}\t{:.2}\t{:.4}", eta, n, alpha, put);
            }
        }
    }
} 