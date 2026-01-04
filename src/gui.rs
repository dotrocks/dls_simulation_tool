use crate::report::export_csv;
use crate::report::export_pdf;
use crate::{
    simulation::simulate,
    structs::{DLSResult, SimulationParams},
};
use chrono::Utc;
use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};
use rfd;
use std::sync::{Arc, Mutex};

pub struct DLSApp {
    params: SimulationParams,
    results: Arc<Mutex<Option<DLSResult>>>,
    loading: bool,
}

impl Default for DLSApp {
    fn default() -> Self {
        Self {
            params: SimulationParams::default(),
            results: Arc::new(Mutex::new(None)),
            loading: false,
        }
    }
}

impl DLSApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self::default()
    }

    fn run_simulation(&mut self, ctx: &egui::Context) {
        self.loading = true;
        let params = self.params.clone();
        let ctx = ctx.clone();
        let results = Arc::clone(&self.results);

        std::thread::spawn(move || {
            let result = simulate(params);
            let mut results_lock = results.lock().unwrap();
            *results_lock = Some(result);
            ctx.request_repaint();
        });
    }
}

impl eframe::App for DLSApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("header")
            .exact_height(70.0)
            .show(ctx, |ui| {
                ui.add_space(10.0);
                ui.vertical_centered(|ui| {
                    ui.label(
                        egui::RichText::new("DYNAMIC LIGHT SCATTERING")
                            .size(26.0)
                            .strong()
                            .color(egui::Color32::from_rgb(52, 152, 219)),
                    );
                });
            });

        egui::SidePanel::left("params")
            .min_width(260.0)
            .max_width(260.0)
            .show(ctx, |ui| {
                ui.add_space(3.0);
                ui.label(
                    egui::RichText::new("PARAMETERS")
                        .size(13.0)
                        .strong()
                        .color(egui::Color32::from_rgb(230, 126, 34)),
                );
                ui.separator();

                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.add(
                        egui::Slider::new(&mut self.params.total_time, 0.001..=1.0)
                            .text("Time (s)")
                            .custom_formatter(|n, _| format!("{:.3}", n)),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.dt, 1e-7..=1e-3)
                            .logarithmic(true)
                            .text("dt (s)")
                            .custom_formatter(|n, _| format!("{:.1e}", n)),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.mean_size_nm, 10.0..=500.0)
                            .text("Mean (nm)")
                            .custom_formatter(|n, _| format!("{:.1}", n)),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.std_size_nm, 1.0..=100.0)
                            .text("Std (nm)")
                            .custom_formatter(|n, _| format!("{:.1}", n)),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.n_particles, 10..=1000)
                            .text("Particles"),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.temperature_c, 0.0..=100.0)
                            .text("Temp (°C)")
                            .custom_formatter(|n, _| format!("{:.1}", n)),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.viscosity_mpa_s, 0.1..=10.0)
                            .text("Viscosity")
                            .custom_formatter(|n, _| format!("{:.2}", n)),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.wavelength_nm, 400.0..=800.0)
                            .text("λ (nm)")
                            .custom_formatter(|n, _| format!("{:.0}", n)),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.scattering_angle_deg, 30.0..=173.0)
                            .text("Angle (°)")
                            .custom_formatter(|n, _| format!("{:.0}", n)),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.beta, 0.1..=1.0)
                            .text("Coherence")
                            .custom_formatter(|n, _| format!("{:.2}", n)),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.shot_noise_level, 0.0..=1.0)
                            .text("Shot noise")
                            .custom_formatter(|n, _| format!("{:.2}", n)),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.detector_noise_level, 0.0..=0.5)
                            .text("Detector")
                            .custom_formatter(|n, _| format!("{:.2}", n)),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.dark_count_rate, 0.0..=5000.0)
                            .text("Dark count")
                            .custom_formatter(|n, _| format!("{:.0}", n)),
                    );
                });

                ui.add_space(6.0);
                ui.separator();

                let run_btn = egui::Button::new(egui::RichText::new("RUN").size(13.0).strong())
                    .fill(egui::Color32::from_rgb(46, 204, 113));
                if ui
                    .add_sized([ui.available_width(), 32.0], run_btn)
                    .clicked()
                {
                    self.run_simulation(ctx);
                }

                if ui
                    .add_sized(
                        [ui.available_width(), 24.0],
                        egui::Button::new("Reset").small(),
                    )
                    .clicked()
                {
                    self.params = SimulationParams::default();
                }

                if self.loading {
                    ui.add_space(6.0);
                    ui.horizontal(|ui| {
                        ui.spinner();
                        ui.label(
                            egui::RichText::new("Computing...")
                                .size(11.0)
                                .color(egui::Color32::from_rgb(52, 152, 219)),
                        );
                    });
                }

                ui.add_space(6.0);
                ui.separator();
                ui.label(
                    egui::RichText::new("RESULTS")
                        .size(13.0)
                        .strong()
                        .color(egui::Color32::from_rgb(155, 89, 182)),
                );

                let results_lock = self.results.lock().unwrap();
                if let Some(results) = results_lock.as_ref() {
                    self.loading = false;

                    ui.add_space(3.0);
                    egui::Grid::new("stats_grid")
                        .num_columns(2)
                        .spacing([6.0, 2.0])
                        .show(ui, |ui| {
                            ui.label(egui::RichText::new("SNR:").size(11.0).strong());
                            ui.label(
                                egui::RichText::new(format!(
                                    "{:.1} dB",
                                    results.noise_metrics.snr_db
                                ))
                                .size(11.0)
                                .color(egui::Color32::from_rgb(46, 204, 113)),
                            );
                            ui.end_row();

                            ui.label(egui::RichText::new("Mean:").size(11.0).strong());
                            ui.label(
                                egui::RichText::new(format!(
                                    "{:.1} nm",
                                    results.stats.mean_size_int
                                ))
                                .size(11.0),
                            );
                            ui.end_row();

                            ui.label(egui::RichText::new("Std:").size(11.0).strong());
                            ui.label(
                                egui::RichText::new(format!(
                                    "{:.1} nm",
                                    results.stats.std_size_int
                                ))
                                .size(11.0),
                            );
                            ui.end_row();

                            ui.label(egui::RichText::new("PDI:").size(11.0).strong());
                            ui.label(
                                egui::RichText::new(format!(
                                    "{:.3}",
                                    results.stats.polydispersity_int
                                ))
                                .size(11.0),
                            );
                            ui.end_row();
                        });

                    ui.add_space(6.0);

                    if ui
                        .add_sized(
                            [ui.available_width(), 24.0],
                            egui::Button::new("Export PDF").small(),
                        )
                        .clicked()
                    {
                        let default_name =
                            format!("dls_report_{}.pdf", Utc::now().format("%Y-%m-%d_%H-%M-%S"));
                        if let Some(path) = rfd::FileDialog::new()
                            .add_filter("PDF Files", &["pdf"])
                            .set_file_name(&default_name)
                            .save_file()
                        {
                            if let Err(err) =
                                export_pdf(&path.to_string_lossy(), &self.params, results)
                            {
                                eprintln!("PDF export error: {err}");
                            }
                        }
                    }

                    if ui
                        .add_sized(
                            [ui.available_width(), 24.0],
                            egui::Button::new("Export CSVs").small(),
                        )
                        .clicked()
                    {
                        if let Some(path) = rfd::FileDialog::new().pick_folder() {
                            if let Err(err) =
                                export_csv(&path.to_string_lossy(), &self.params, results)
                            {
                                eprintln!("CSV export error: {err}");
                            }
                        }
                    }
                } else {
                    ui.add_space(3.0);
                    ui.label(
                        egui::RichText::new("No data yet...")
                            .size(10.0)
                            .italics()
                            .color(egui::Color32::DARK_GRAY),
                    );
                }
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.add_space(8.0);

                let results_lock = self.results.lock().unwrap();
                if let Some(results) = results_lock.as_ref() {
                    let available_width = ui.available_width();
                    let plot_width = available_width - 20.0;
                    let plot_height = 220.0;

                    ui.group(|ui| {
                        ui.label(
                            egui::RichText::new("SCATTERED INTENSITY vs TIME")
                                .size(14.0)
                                .strong()
                                .color(egui::Color32::from_rgb(231, 76, 60)),
                        );
                        Plot::new("intensity_plot")
                            .width(plot_width)
                            .height(plot_height)
                            .allow_zoom(true)
                            .allow_drag(true)
                            .x_axis_label("Time (s)")
                            .y_axis_label("Intensity (a.u.)")
                            .legend(
                                egui_plot::Legend::default().position(egui_plot::Corner::RightTop),
                            )
                            .show(ui, |plot_ui| {
                                let n_show = 1000.min(results.time.len());
                                let ideal_points: PlotPoints = (0..n_show)
                                    .map(|i| [results.time[i], results.intensity_ideal[i]])
                                    .collect();
                                plot_ui.line(
                                    Line::new("Ideal Signal", ideal_points)
                                        .name("Ideal Signal")
                                        .color(egui::Color32::from_rgb(46, 204, 113))
                                        .width(2.5),
                                );

                                let noisy_points: PlotPoints = (0..n_show)
                                    .map(|i| [results.time[i], results.intensity_noisy[i]])
                                    .collect();
                                plot_ui.line(
                                    Line::new("Noisy Signal", noisy_points)
                                        .name("Noisy Signal")
                                        .color(egui::Color32::from_rgb(231, 76, 60))
                                        .width(1.8),
                                );
                            });
                    });

                    ui.add_space(12.0);

                    ui.group(|ui| {
                        ui.label(
                            egui::RichText::new("FIELD AUTOCORRELATION FUNCTION g₁(τ)")
                                .size(14.0)
                                .strong()
                                .color(egui::Color32::from_rgb(52, 152, 219)),
                        );
                        Plot::new("g1_plot")
                            .width(plot_width)
                            .height(plot_height)
                            .allow_zoom(true)
                            .allow_drag(true)
                            .x_axis_label("log₁₀(τ) [s]")
                            .y_axis_label("g₁(τ)")
                            .legend(
                                egui_plot::Legend::default().position(egui_plot::Corner::RightTop),
                            )
                            .show(ui, |plot_ui| {
                                let g1_sim: PlotPoints = results
                                    .tau
                                    .iter()
                                    .zip(&results.g1_numeric_ideal)
                                    .filter(|(tau, _)| **tau > 1e-10)
                                    .map(|(tau, g1)| [tau.log10(), *g1])
                                    .collect();
                                plot_ui.line(
                                    Line::new("Simulated", g1_sim)
                                        .name("Simulated")
                                        .color(egui::Color32::from_rgb(231, 76, 60))
                                        .width(2.5),
                                );

                                let g1_theory: PlotPoints = results
                                    .tau_theory
                                    .iter()
                                    .zip(&results.g1_theory)
                                    .filter(|(tau, _)| **tau > 1e-10)
                                    .map(|(tau, g1)| [tau.log10(), *g1])
                                    .collect();
                                plot_ui.line(
                                    Line::new("Theoretical", g1_theory)
                                        .name("Theoretical")
                                        .color(egui::Color32::from_rgb(46, 204, 113))
                                        .width(2.5),
                                );
                            });
                    });

                    ui.add_space(12.0);

                    ui.group(|ui| {
                        ui.label(
                            egui::RichText::new("INTENSITY AUTOCORRELATION FUNCTION g₂(τ)")
                                .size(14.0)
                                .strong()
                                .color(egui::Color32::from_rgb(155, 89, 182)),
                        );
                        Plot::new("g2_plot")
                            .width(plot_width)
                            .height(plot_height)
                            .allow_zoom(true)
                            .allow_drag(true)
                            .x_axis_label("log₁₀(τ) [s]")
                            .y_axis_label("g₂(τ)")
                            .legend(
                                egui_plot::Legend::default().position(egui_plot::Corner::RightTop),
                            )
                            .show(ui, |plot_ui| {
                                let g2_sim: PlotPoints = results
                                    .tau
                                    .iter()
                                    .zip(&results.g2_numeric_noisy)
                                    .filter(|(tau, _)| **tau > 1e-10)
                                    .map(|(tau, g2)| [tau.log10(), *g2])
                                    .collect();
                                plot_ui.line(
                                    Line::new("Simulated (Noisy)", g2_sim)
                                        .name("Simulated (Noisy)")
                                        .color(egui::Color32::from_rgb(231, 76, 60))
                                        .width(2.5),
                                );

                                let g2_theory: PlotPoints = results
                                    .tau_theory
                                    .iter()
                                    .zip(&results.g2_theory)
                                    .filter(|(tau, _)| **tau > 1e-10)
                                    .map(|(tau, g2)| [tau.log10(), *g2])
                                    .collect();
                                plot_ui.line(
                                    Line::new("Theoretical", g2_theory)
                                        .name("Theoretical")
                                        .color(egui::Color32::from_rgb(46, 204, 113))
                                        .width(2.5),
                                );
                            });
                    });

                    ui.add_space(12.0);

                    ui.group(|ui| {
                        ui.label(
                            egui::RichText::new("PARTICLE SIZE DISTRIBUTION")
                                .size(14.0)
                                .strong()
                                .color(egui::Color32::from_rgb(230, 126, 34)),
                        );
                        Plot::new("size_plot")
                            .width(plot_width)
                            .height(plot_height)
                            .allow_zoom(true)
                            .allow_drag(true)
                            .x_axis_label("Size (nm)")
                            .y_axis_label("Relative Intensity")
                            .legend(
                                egui_plot::Legend::default().position(egui_plot::Corner::RightTop),
                            )
                            .show(ui, |plot_ui| {
                                let int_dist: PlotPoints = results
                                    .size_nm
                                    .iter()
                                    .zip(&results.size_intensity)
                                    .map(|(size, intensity)| [*size, *intensity])
                                    .collect();
                                plot_ui.line(
                                    Line::new("Intensity Weighted", int_dist)
                                        .name("Intensity Weighted")
                                        .color(egui::Color32::from_rgb(155, 89, 182))
                                        .width(3.0),
                                );

                                let num_dist: PlotPoints = results
                                    .size_num_nm
                                    .iter()
                                    .zip(&results.size_num_dist)
                                    .map(|(size, num)| [*size, *num])
                                    .collect();
                                plot_ui.line(
                                    Line::new("Number Based", num_dist)
                                        .name("Number Based")
                                        .color(egui::Color32::from_rgb(243, 156, 18))
                                        .width(3.0),
                                );
                            });
                    });

                    ui.add_space(15.0);
                } else {
                    ui.vertical_centered(|ui| {
                        ui.add_space(120.0);
                        ui.label(egui::RichText::new("🔬").size(64.0));
                        ui.add_space(15.0);
                        ui.label(
                            egui::RichText::new("AWAITING SIMULATION DATA")
                                .size(18.0)
                                .strong()
                                .color(egui::Color32::GRAY),
                        );
                        ui.add_space(8.0);
                        ui.label(
                            egui::RichText::new("Configure parameters and click 'RUN SIMULATION'")
                                .size(13.0)
                                .color(egui::Color32::DARK_GRAY),
                        );
                    });
                }
            });
        });
    }
}
