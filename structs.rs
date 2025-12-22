use serde::{Deserialize, Serialize};

#[derive(Clone, Serialize, Deserialize)]
pub struct SimulationParams {
    /// Toplam ölçüm süresi (s)
    pub total_time: f64,
    /// Nominal zaman adımı (s) – simülasyonda üst limit ile clamp ediliyor
    pub dt: f64,
    /// Boyut dağılımının ortalaması (nm) – intensity-weighted lognormal mean
    pub mean_size_nm: f64,
    /// Boyut dağılımının standart sapması (nm)
    pub std_size_nm: f64,
    /// Simüle edilen parçacık sayısı
    pub n_particles: usize,
    /// Sıcaklık (°C)
    pub temperature_c: f64,
    /// Ortam viskozitesi (mPa·s)
    pub viscosity_mpa_s: f64,
    /// Lazer dalga boyu (nm)
    pub wavelength_nm: f64,
    /// Saçılma açısı (derece)
    pub scattering_angle_deg: f64,
    /// Koherens faktörü (β) – Siegert ilişkisi için
    pub beta: f64,
    /// Shot noise ölçek faktörü (0–1)
    pub shot_noise_level: f64,
    /// Dedektör gürültü seviyesi (Gauss, birim: I_mean)
    pub detector_noise_level: f64,
    /// Dark count oranı (counts/s)
    pub dark_count_rate: f64,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct DLSResult {
    /// Zaman ekseni (s)
    pub time: Vec<f64>,
    /// Gürültüsüz saçılan yoğunluk (I_ideal)
    pub intensity_ideal: Vec<f64>,
    /// Gürültü eklenmiş yoğunluk (I_noisy)
    pub intensity_noisy: Vec<f64>,

    /// Sayısal g1 ve g2 için kullanılan lag zamanları (s)
    pub tau: Vec<f64>,
    /// Sayısal olarak hesaplanmış field otokorelasyonu g1(τ)
    pub g1_numeric_ideal: Vec<f64>,
    /// g1'den Siegert ile türetilmiş intensity otokorelasyonu g2(τ)
    pub g2_numeric_noisy: Vec<f64>,

    /// Teorik g1 ve g2 için lag zamanları (s)
    pub tau_theory: Vec<f64>,
    /// Teorik g1(τ) (çoklu boyutlu, ağırlıklandırılmış exponential)
    pub g1_theory: Vec<f64>,
    /// Teorik g2(τ) = 1 + β g1(τ)²
    pub g2_theory: Vec<f64>,

    /// Intensity-weighted boyut dağılımı için x ekseni (nm)
    pub size_nm: Vec<f64>,
    /// Intensity-weighted boyut dağılımı (normalize)
    pub size_intensity: Vec<f64>,

    /// Number-based boyut dağılımı için x ekseni (nm)
    pub size_num_nm: Vec<f64>,
    /// Number-based boyut dağılımı (normalize)
    pub size_num_dist: Vec<f64>,

    /// Tek tek çekilmiş ham partikül boyutları (nm)
    pub raw_sizes_nm: Vec<f64>,
    /// Her partikül için intensity ağırlığı (∝ r^6, normalize)
    pub intensity_weights: Vec<f64>,

    pub noise_metrics: NoiseMetrics,
    pub stats: Stats,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct NoiseMetrics {
    /// Sinyal/gürültü oranı (dB)
    pub snr_db: f64,
    /// Gürültü gücü (variance of noisy - ideal)
    pub noise_power: f64,
    /// Sinyal gücü (variance of I_ideal)
    pub signal_power: f64,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Stats {
    /// Intensity-weighted ortalama boyut (Z-average'e yakın anlam) [nm]
    pub mean_size_int: f64,
    pub std_size_int: f64,
    /// Number-based ortalama boyut [nm]
    pub mean_size_num: f64,
    pub std_size_num: f64,
    /// Intensity-weighted PDI (σ/μ) – 0–1 aralığına yakın beklenir [web:3][web:86]
    pub polydispersity_int: f64,
    /// Number-based PDI (σ/μ)
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
