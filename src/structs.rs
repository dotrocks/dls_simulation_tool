use serde::{Deserialize, Serialize};

#[derive(Clone, Serialize, Deserialize)]
pub struct SimulationParams {
    pub total_time: f64,
    pub dt: f64,
    pub mean_size_nm: f64,
    pub std_size_nm: f64,
    pub n_particles: usize,
    pub temperature_c: f64,
    pub viscosity_mpa_s: f64,
    pub wavelength_nm: f64,
    pub scattering_angle_deg: f64,
    pub beta: f64,
    pub shot_noise_level: f64,
    pub detector_noise_level: f64,
    pub dark_count_rate: f64,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct DLSResult {
    pub time: Vec<f64>,
    pub intensity_ideal: Vec<f64>,
    pub intensity_noisy: Vec<f64>,
    pub tau: Vec<f64>,
    pub g1_numeric_ideal: Vec<f64>,
    pub g2_numeric_noisy: Vec<f64>,
    pub tau_theory: Vec<f64>,
    pub g1_theory: Vec<f64>,
    pub g2_theory: Vec<f64>,
    pub size_nm: Vec<f64>,
    pub size_intensity: Vec<f64>,
    pub size_num_nm: Vec<f64>,
    pub size_num_dist: Vec<f64>,
    pub raw_sizes_nm: Vec<f64>,
    pub intensity_weights: Vec<f64>,
    pub noise_metrics: NoiseMetrics,
    pub stats: Stats,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct NoiseMetrics {
    pub snr_db: f64,
    pub noise_power: f64,
    pub signal_power: f64,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Stats {
    pub mean_size_int: f64,
    pub std_size_int: f64,
    pub mean_size_num: f64,
    pub std_size_num: f64,
    pub polydispersity_int: f64,
    pub polydispersity_num: f64,
}

impl Default for SimulationParams {
    fn default() -> Self {
        Self {
            total_time: 0.1,
            dt: 1e-5,
            mean_size_nm: 100.0,
            std_size_nm: 20.0,
            n_particles: 200,
            temperature_c: 25.0,
            viscosity_mpa_s: 0.89,
            wavelength_nm: 633.0,
            scattering_angle_deg: 90.0,
            beta: 0.8,
            shot_noise_level: 0.1,
            detector_noise_level: 0.05,
            dark_count_rate: 1000.0,
        }
    }
}
