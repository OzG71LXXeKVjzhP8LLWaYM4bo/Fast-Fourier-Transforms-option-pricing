# Fast-Fourier-Transforms-option-pricing

This repository provides Fast Fourier Transform (FFT) pricing for European options using the Carr–Madan method with three models:
- Black–Scholes
- Heston (stochastic volatility)
- Variance–Gamma

## Crate layout
- `fft/`: Rust crate with a library (`src/lib.rs`) and three binaries in `src/bin/`:
  - `bs.rs` (Black–Scholes)
  - `heston.rs` (Heston)
  - `vg.rs` (Variance–Gamma)

## Build
```bash
cd fft
cargo build
```

## Run
Example (defaults match the article’s style):
```bash
cargo run --bin bs
cargo run --bin heston
cargo run --bin vg
```
Optional args (shared): `--s0`, `--k`, `--r`, `--q`, `--t`.
Model-specific:
- BS: `--sigma`
- Heston: `--kappa`, `--theta`, `--vol-of-vol`, `--rho`, `--v0`
- VG: `--sigma`, `--nu`, `--theta`

Each binary prints a parameter sweep table with columns: eta, 2^n, alpha, put.

## Reference
This implementation follows the Carr–Madan FFT approach and mirrors the walkthrough from the Medium article “Pricing Options with Fast Fourier Transform (FFT) — Heston, Variance-Gamma and Black-Scholes (with code)” by Antoni Smolski: [Medium article](https://antonismolski.medium.com/pricing-options-with-fast-fourier-transform-fft-ee7ae0a97c43)