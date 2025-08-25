use clap::Parser;
use fft::fft_pricer::{CarrMadanParams, price_put_at_strike};
use fft::models::Heston;

#[derive(Parser, Debug)]
#[command(name = "heston", about = "Heston put pricing via Carr–Madan FFT")] 
struct Args {
    #[arg(long, default_value = "100.0")] s0: f64,
    #[arg(long, default_value = "80.0")] k: f64,
    #[arg(long, default_value = "0.055")] r: f64,
    #[arg(long, default_value = "0.03")] q: f64,
    #[arg(long, default_value = "1.0")] t: f64,
    #[arg(long, default_value = "2.0")] kappa: f64,
    #[arg(long, default_value = "0.05")] theta: f64,
    #[arg(long, default_value = "0.3")] vol_of_vol: f64,
    #[arg(long, default_value = "-0.7")] rho: f64,
    #[arg(long, default_value = "0.04")] v0: f64,
}

fn main() {
    let args = Args::parse();
    let model = Heston { kappa: args.kappa, theta: args.theta, vol_of_vol: args.vol_of_vol, rho: args.rho, v0: args.v0 };

    let alpha_values = vec![1.01, 1.25, 1.50, 1.75, 2.00, 5.00];
    let eta_values = vec![0.10, 0.25];
    let n_values = vec![6usize, 10usize];

    println!("Model = Heston");
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