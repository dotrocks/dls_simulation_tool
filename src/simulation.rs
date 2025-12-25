use crate::structs::{DLSResult, NoiseMetrics, SimulationParams, Stats};
use rand::{Rng, rng, rngs::ThreadRng};
use rand_distr::{Distribution, LogNormal, Normal};
use std::f64::consts::PI;

const K_B: f64 = 1.380649e-23;

pub struct DLSSimulator {
    rng: ThreadRng,
    normal_dist: Normal<f64>,
}

impl DLSSimulator {
    pub fn new() -> Self {
        Self {
            rng: rng(),
            normal_dist: Normal::new(0.0, 1.0).unwrap(),
        }
    }

    pub fn simulate_dls(&mut self, params: SimulationParams) -> DLSResult {
    let mut n_steps = (params.total_time / params.dt).ceil() as usize;
    let max_steps = 1_000_000;
    if n_steps > max_steps {
        n_steps = max_steps;
    }
    let effective_dt = params.total_time / n_steps as f64;
    let time: Vec<f64> = (0..n_steps).map(|i| (i as f64) * effective_dt).collect();


        let t = params.temperature_c + 273.15;
        let eta = params.viscosity_mpa_s * 1e-3;
        let wavelength_m = params.wavelength_nm * 1e-9;
        let n_medium = 1.33;
        let theta = params.scattering_angle_deg.to_radians();
        let q = (4.0 * PI * n_medium / wavelength_m) * (theta / 2.0).sin();
        let cv = (params.std_size_nm / params.mean_size_nm).max(0.01);
        let log_std = (1.0 + cv.powi(2)).ln().sqrt();
        let log_mean = (params.mean_size_nm / (1.0 + cv.powi(2)).sqrt()).ln();


        let log_normal = LogNormal::new(log_mean, log_std).unwrap();
        let sizes_nm: Vec<f64> = (0..params.n_particles)
            .map(|_| {
                let sample = log_normal.sample(&mut self.rng);
                sample.max(1.0)
            })
            .collect();

        let radii_m: Vec<f64> = sizes_nm.iter().map(|s| (s / 2.0) * 1e-9).collect();
        let d: Vec<f64> = radii_m
            .iter()
            .map(|r| K_B * t / (6.0 * PI * eta * r))
            .collect();
        let gamma: Vec<f64> = d.iter().map(|d_val| d_val * q.powi(2)).collect();

        let mut intensity_weights: Vec<f64> = radii_m.iter().map(|r| r.powi(6)).collect();
        let sum_weights: f64 = intensity_weights.iter().sum();
        for weight in intensity_weights.iter_mut() {
            *weight = *weight / sum_weights;
        }

        let mut e_total = vec![0.0; n_steps];
        for i in 0..params.n_particles {
            let random_value: f64 = self.rng.random();
            let initial_phase = 2.0 * PI * random_value;
            let sqrt_2_gamma = (2.0 * gamma[i]).sqrt();
            let mut phi = initial_phase;

            for j in 0..n_steps {
                let normal_rv = self.normal_dist.sample(&mut self.rng);
                let dw = normal_rv * effective_dt.sqrt();
                phi = phi + sqrt_2_gamma * dw;
                e_total[j] = e_total[j] + intensity_weights[i] * phi.cos();
            }
        }

        let i_ideal: Vec<f64> = e_total.iter().map(|e| e.powi(2)).collect();
        let i_mean = i_ideal.iter().sum::<f64>() / (i_ideal.len() as f64);

        let mut i_final = vec![0.0; n_steps];
        for (j, t_val) in time.iter().enumerate() {
            let photon_count = i_ideal[j] / i_mean * 1000.0;
            let shot_noise = self.sample_poisson(photon_count * params.shot_noise_level);
            let mut intensity = i_ideal[j]
                + (shot_noise - photon_count * params.shot_noise_level) * i_mean / 1000.0;

            let detector_noise =
                self.normal_dist.sample(&mut self.rng) * params.detector_noise_level * i_mean;
            intensity = intensity + detector_noise;

            let dark_noise = self.sample_poisson(params.dark_count_rate * effective_dt);
            intensity += (dark_noise - params.dark_count_rate * effective_dt) * i_mean / 1000.0;


            let time_hours = t_val / 3600.0;
            let drift_coeff = 0.001;
            let baseline_drift = drift_coeff * i_mean * time_hours;
            i_final[j] = (intensity + baseline_drift).max(0.0);

        }

        let noise_signal: Vec<f64> = i_final
            .iter()
            .zip(&i_ideal)
            .map(|(noisy, ideal)| noisy - ideal)
            .collect();
        let noise_power = Self::variance(&noise_signal);
        let signal_power = Self::variance(&i_ideal);
        let snr_db = if noise_power > 0.0 {
            10.0 * (signal_power / noise_power).log10()
        } else {
            100.0
        };

        let max_lag = (n_steps / 4).min(2000);
        let e_mean = e_total.iter().sum::<f64>() / (e_total.len() as f64);
        let e_fluct: Vec<f64> = e_total.iter().map(|e| e - e_mean).collect();

        // Field fluctuation → g1
        let g1_numeric_ideal = Self::compute_autocorr(&e_fluct, max_lag);

        // Siegert relation: g2 = 1 + β |g1|²
        let g2_numeric_noisy: Vec<f64> = g1_numeric_ideal
            .iter()
            .map(|g1| 1.0 + params.beta * g1.powi(2))
            .collect();


        let tau: Vec<f64> = (0..max_lag).map(|i| (i as f64) * effective_dt).collect();


        let mut g1_theory = vec![0.0; max_lag];
        for (lag_idx, tau_val) in tau.iter().enumerate() {
            let mut sum_val = 0.0;
            for i in 0..params.n_particles {
                sum_val = sum_val + intensity_weights[i] * (-gamma[i] * tau_val).exp();
            }
            g1_theory[lag_idx] = sum_val;
        }
        let g2_theory: Vec<f64> = g1_theory
            .iter()
            .map(|g1| 1.0 + params.beta * g1.powi(2))
            .collect();

        let bins = 30.min(params.n_particles / 10).max(1);
        let (hist_intensity, hist_edges) =
            Self::weighted_histogram(&sizes_nm, &intensity_weights, bins);
        let hist_x: Vec<f64> = hist_edges.windows(2).map(|w| (w[0] + w[1]) / 2.0).collect();
        let (hist_num, edges_num) = Self::simple_histogram(&sizes_nm, bins);
        let hist_num_x: Vec<f64> = edges_num.windows(2).map(|w| (w[0] + w[1]) / 2.0).collect();

        let (mean_size_int, std_size_int) = Self::weighted_stats(&sizes_nm, &intensity_weights);
        let mean_size_num = sizes_nm.iter().sum::<f64>() / (sizes_nm.len() as f64);
        let variance_num = sizes_nm
            .iter()
            .map(|s| (s - mean_size_num).powi(2))
            .sum::<f64>()
            / (sizes_nm.len() as f64);
        let std_size_num = variance_num.sqrt();

        DLSResult {
            time,
            intensity_ideal: i_ideal,
            intensity_noisy: i_final,
            tau: tau.clone(),
            g1_numeric_ideal,
            g2_numeric_noisy,
            tau_theory: tau,
            g1_theory,
            g2_theory,
            size_nm: hist_x,
            size_intensity: hist_intensity,
            size_num_nm: hist_num_x,
            size_num_dist: hist_num,
            raw_sizes_nm: sizes_nm,
            intensity_weights,
            noise_metrics: NoiseMetrics {
                snr_db,
                noise_power,
                signal_power,
            },
            stats: Stats {
                mean_size_int,
                std_size_int,
                mean_size_num,
                std_size_num,
                polydispersity_int: std_size_int / mean_size_int,
                polydispersity_num: std_size_num / mean_size_num,
            },
        }
    }

    fn sample_poisson(&mut self, lambda: f64) -> f64 {
        if lambda < 30.0 {
            let l = (-lambda).exp();
            let mut k = 0;
            let mut p = 1.0;
            while p > l {
                k = k + 1;
                let random_val: f64 = self.rng.random();
                p = p * random_val;
            }
            (k - 1) as f64
        } else {
            let normal_dist = Normal::new(lambda, lambda.sqrt()).unwrap();
            let sample = normal_dist.sample(&mut self.rng);
            sample.max(0.0)
        }
    }

    fn variance(data: &[f64]) -> f64 {
        let mean = data.iter().sum::<f64>() / (data.len() as f64);
        data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (data.len() as f64)
    }

    fn compute_autocorr(data: &[f64], max_lag: usize) -> Vec<f64> {
        let n = data.len();
        let mean = data.iter().sum::<f64>() / n as f64;
        let centered: Vec<f64> = data.iter().map(|x| x - mean).collect();

        let var = centered
            .iter()
            .map(|x| x * x)
            .sum::<f64>() / n as f64;

        let mut result = vec![0.0; max_lag];
        for lag in 0..max_lag.min(n) {
            let mut sum = 0.0;
            for i in 0..(n - lag) {
                sum += centered[i] * centered[i + lag];
            }
            // normalize by N and variance → g(0)=1
            result[lag] = sum / (n as f64 * var);
        }
        result
    }


    fn weighted_stats(data: &[f64], weights: &[f64]) -> (f64, f64) {
        let sum_w: f64 = weights.iter().sum();
        let mean = data.iter().zip(weights).map(|(x, w)| x * w).sum::<f64>() / sum_w;
        let variance = data
            .iter()
            .zip(weights)
            .map(|(x, w)| w * (x - mean).powi(2))
            .sum::<f64>()
            / sum_w;
        (mean, variance.sqrt())
    }

    fn weighted_histogram(data: &[f64], weights: &[f64], bins: usize) -> (Vec<f64>, Vec<f64>) {
        if data.is_empty() {
            return (vec![0.0; bins], vec![0.0; bins + 1]);
        }
        let min = data.iter().fold(f64::INFINITY, |a, b| a.min(*b));
        let max = data.iter().fold(f64::NEG_INFINITY, |a, b| a.max(*b));
        let bin_width = (max - min) / (bins as f64);
        let mut histogram = vec![0.0; bins];
        let edges: Vec<f64> = (0..=bins).map(|i| min + (i as f64) * bin_width).collect();

        for (x, w) in data.iter().zip(weights) {
            let bin_idx = ((x - min) / bin_width).floor() as usize;
            if bin_idx < bins {
                histogram[bin_idx] = histogram[bin_idx] + w;
            }
        }
        let sum: f64 = histogram.iter().sum();
        if sum > 0.0 {
            for h in histogram.iter_mut() {
                *h = *h / sum;
            }
        }
        (histogram, edges)
    }

    fn simple_histogram(data: &[f64], bins: usize) -> (Vec<f64>, Vec<f64>) {
        if data.is_empty() {
            return (vec![0.0; bins], vec![0.0; bins + 1]);
        }
        let min = data.iter().fold(f64::INFINITY, |a, b| a.min(*b));
        let max = data.iter().fold(f64::NEG_INFINITY, |a, b| a.max(*b));
        let bin_width = (max - min) / (bins as f64);
        let mut histogram = vec![0.0; bins];
        let edges: Vec<f64> = (0..=bins).map(|i| min + (i as f64) * bin_width).collect();

        for x in data {
            let bin_idx = ((x - min) / bin_width).floor() as usize;
            if bin_idx < bins {
                histogram[bin_idx] = histogram[bin_idx] + 1.0;
            }
        }
        let sum = data.len() as f64;
        for h in histogram.iter_mut() {
            *h = *h / sum;
        }
        (histogram, edges)
    }
}

pub fn simulate(params: SimulationParams) -> DLSResult {
    let mut simulator = DLSSimulator::new();
    simulator.simulate_dls(params)
}
