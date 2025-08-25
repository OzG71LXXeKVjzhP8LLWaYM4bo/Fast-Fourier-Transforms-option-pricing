use clap::Parser;
use fft::fft_pricer::{CarrMadanParams, price_put_at_strike};
use fft::models::BlackScholes;

#[derive(Parser, Debug)]
#[command(name = "bs", about = "Black-Scholes put pricing via Carr–Madan FFT")] 
struct Args {
    #[arg(long, default_value = "100.0")] s0: f64,
    #[arg(long, default_value = "80.0")] k: f64,
    #[arg(long, default_value = "0.055")] r: f64,
    #[arg(long, default_value = "0.03")] q: f64,
    #[arg(long, default_value = "1.0")] t: f64,
    #[arg(long, default_value = "0.3")] sigma: f64,
    /// Optional FFT damping (if omitted, print sweep)
    #[arg(long)] alpha: Option<f64>,
    /// Optional FFT step size (if omitted, print sweep)
    #[arg(long)] eta: Option<f64>,
    /// Optional FFT exponent n, N = 2^n (if omitted, print sweep)
    #[arg(long)] n: Option<usize>,
}

fn main() {
    let args = Args::parse();
    let model = BlackScholes { sigma: args.sigma };

    let alpha_values = match args.alpha { Some(a) => vec![a], None => vec![1.01, 1.25, 1.50, 1.75, 2.00, 5.00] };
    let eta_values = match args.eta { Some(e) => vec![e], None => vec![0.10, 0.25] };
    let n_values = match args.n { Some(n) => vec![n], None => vec![6usize, 10usize] };

    println!("Model = BS");
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