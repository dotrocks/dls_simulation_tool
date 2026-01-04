use crate::structs::{DLSResult, SimulationParams};
use chrono::Utc;
use csv;
use genpdf::Element;
use genpdf::{Alignment, elements, style};
use opener;
use plotters::prelude::*;
use std::fs::File;
use std::path::Path;
use tempfile::NamedTempFile;

macro_rules! row {
    ($table:ident, $label:expr, $value:expr) => {{
        $table
            .row()
            .element(elements::Paragraph::new($label))
            .element(elements::Paragraph::new($value))
            .push()
            .ok();
    }};
}

pub fn export_pdf(path: &str, params: &SimulationParams, result: &DLSResult) -> anyhow::Result<()> {
    let fonts_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("font")
        .join("fonts")
        .join("Frozen Fonts")
        .join("Monaspace Xenon");

    let font_family = genpdf::fonts::from_files(&fonts_dir, "MonaspaceXenonFrozen", None)
        .expect("Failed to load font family");

    let mut doc = genpdf::Document::new(font_family);
    doc.set_title("DLS Simulation Report");

    let mut decorator = genpdf::SimplePageDecorator::new();
    decorator.set_margins(20);
    doc.set_page_decorator(decorator);

    let heading = elements::Paragraph::new("Dynamic Light Scattering Simulation")
        .aligned(Alignment::Center)
        .styled(style::Style::new().bold().with_font_size(18));
    doc.push(heading);
    doc.push(elements::Break::new(1));

    doc.push(
        elements::Paragraph::new("Simulation Parameters")
            .styled(style::Style::new().bold().with_font_size(14)),
    );
    doc.push(elements::Break::new(0.5));

    let mut table = elements::TableLayout::new(vec![1, 1]);

    row!(table, "Total time (s)", format!("{:.4}", params.total_time));
    row!(table, "dt (s)", format!("{:.1e}", params.dt));
    row!(
        table,
        "Mean size (nm)",
        format!("{:.1}", params.mean_size_nm)
    );
    row!(table, "Std size (nm)", format!("{:.1}", params.std_size_nm));
    row!(table, "Particles", format!("{}", params.n_particles));
    row!(
        table,
        "Temperature (°C)",
        format!("{:.1}", params.temperature_c)
    );
    row!(
        table,
        "Viscosity (mPa·s)",
        format!("{:.2}", params.viscosity_mpa_s)
    );
    row!(
        table,
        "Wavelength (nm)",
        format!("{:.0}", params.wavelength_nm)
    );
    row!(
        table,
        "Angle (°)",
        format!("{:.0}", params.scattering_angle_deg)
    );
    row!(table, "Beta", format!("{:.2}", params.beta));
    row!(
        table,
        "Shot noise",
        format!("{:.2}", params.shot_noise_level)
    );
    row!(
        table,
        "Detector noise",
        format!("{:.2}", params.detector_noise_level)
    );
    row!(
        table,
        "Dark count (cps)",
        format!("{:.0}", params.dark_count_rate)
    );

    doc.push(table);
    doc.push(elements::Break::new(1));

    doc.push(
        elements::Paragraph::new("Results Summary")
            .styled(style::Style::new().bold().with_font_size(14)),
    );
    doc.push(elements::Break::new(0.5));

    let mut table2 = elements::TableLayout::new(vec![1, 1]);

    row!(
        table2,
        "SNR (dB)",
        format!("{:.1}", result.noise_metrics.snr_db)
    );
    row!(
        table2,
        "Mean size (intensity, nm)",
        format!("{:.1}", result.stats.mean_size_int)
    );
    row!(
        table2,
        "Std size (intensity, nm)",
        format!("{:.1}", result.stats.std_size_int)
    );
    row!(
        table2,
        "PDI (intensity)",
        format!("{:.3}", result.stats.polydispersity_int)
    );
    row!(
        table2,
        "Mean size (number, nm)",
        format!("{:.1}", result.stats.mean_size_num)
    );
    row!(
        table2,
        "PDI (number)",
        format!("{:.3}", result.stats.polydispersity_num)
    );

    doc.push(table2);

    doc.push(elements::PageBreak::new());

    doc.push(
        elements::Paragraph::new("Plots").styled(style::Style::new().bold().with_font_size(14)),
    );
    doc.push(elements::Break::new(0.5));

    if let Ok(temp_file) = NamedTempFile::new() {
        let image_path = temp_file.path().with_extension("png");
        if generate_g1_plot(&image_path, result).is_ok() {
            if let Ok(img) = elements::Image::from_path(&image_path) {
                doc.push(img.with_alignment(Alignment::Center));
            }
        }
    }
    doc.push(elements::Break::new(0.5));

    if let Ok(temp_file) = NamedTempFile::new() {
        let image_path = temp_file.path().with_extension("png");
        if generate_g2_plot(&image_path, result).is_ok() {
            if let Ok(img) = elements::Image::from_path(&image_path) {
                doc.push(img.with_alignment(Alignment::Center));
            }
        }
    }
    doc.push(elements::Break::new(0.5));

    if let Ok(temp_file) = NamedTempFile::new() {
        let image_path = temp_file.path().with_extension("png");
        if generate_size_plot(&image_path, result).is_ok() {
            if let Ok(img) = elements::Image::from_path(&image_path) {
                doc.push(img.with_alignment(Alignment::Center));
            }
        }
    }

    let file = File::create(path)?;
    doc.render(file)?;

    opener::open(path)?;

    Ok(())
}

fn generate_g1_plot(path: &Path, result: &DLSResult) -> anyhow::Result<()> {
    let root = BitMapBackend::new(path, (2000, 1000)).into_drawing_area();
    root.fill(&WHITE)?;

    let data_sim: Vec<(f64, f64)> = result
        .tau
        .iter()
        .zip(&result.g1_numeric_ideal)
        .filter(|(tau, _)| **tau > 1e-10)
        .map(|(tau, g1)| (tau.log10(), *g1))
        .collect();
    let data_theory: Vec<(f64, f64)> = result
        .tau_theory
        .iter()
        .zip(&result.g1_theory)
        .filter(|(tau, _)| **tau > 1e-10)
        .map(|(tau, g1)| (tau.log10(), *g1))
        .collect();

    let min_x = data_sim
        .iter()
        .map(|&(x, _)| x)
        .fold(f64::INFINITY, f64::min)
        .min(
            data_theory
                .iter()
                .map(|&(x, _)| x)
                .fold(f64::INFINITY, f64::min),
        );
    let max_x = data_sim
        .iter()
        .map(|&(x, _)| x)
        .fold(f64::NEG_INFINITY, f64::max)
        .max(
            data_theory
                .iter()
                .map(|&(x, _)| x)
                .fold(f64::NEG_INFINITY, f64::max),
        );
    let min_y = data_sim
        .iter()
        .chain(&data_theory)
        .map(|&(_, y)| y)
        .fold(f64::INFINITY, f64::min);
    let max_y = data_sim
        .iter()
        .chain(&data_theory)
        .map(|&(_, y)| y)
        .fold(f64::NEG_INFINITY, f64::max);

    let mut chart = ChartBuilder::on(&root)
        .caption("Field Autocorrelation Function g_1(T)", ("sans-serif", 64))
        .margin(10)
        .x_label_area_size(80)
        .y_label_area_size(100)
        .build_cartesian_2d(min_x..max_x, min_y..max_y)?;

    chart
        .configure_mesh()
        .x_desc("log10(T) [s]")
        .y_desc("g_1(T)")
        .label_style(("sans-serif", 32))
        .draw()?;

    chart
        .draw_series(LineSeries::new(data_sim, RED.stroke_width(4)))?
        .label("Simulated")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    chart
        .draw_series(LineSeries::new(data_theory, GREEN.stroke_width(4)))?
        .label("Theoretical")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart
        .configure_series_labels()
        .background_style(&WHITE)
        .label_font(("sans-serif", 32))
        .draw()?;

    root.present()?;
    Ok(())
}

fn generate_g2_plot(path: &Path, result: &DLSResult) -> anyhow::Result<()> {
    let root = BitMapBackend::new(path, (2000, 1000)).into_drawing_area();
    root.fill(&WHITE)?;

    let data_sim: Vec<(f64, f64)> = result
        .tau
        .iter()
        .zip(&result.g2_numeric_noisy)
        .filter(|(tau, _)| **tau > 1e-10)
        .map(|(tau, g2)| (tau.log10(), *g2))
        .collect();
    let data_theory: Vec<(f64, f64)> = result
        .tau_theory
        .iter()
        .zip(&result.g2_theory)
        .filter(|(tau, _)| **tau > 1e-10)
        .map(|(tau, g2)| (tau.log10(), *g2))
        .collect();

    let min_x = data_sim
        .iter()
        .map(|&(x, _)| x)
        .fold(f64::INFINITY, f64::min)
        .min(
            data_theory
                .iter()
                .map(|&(x, _)| x)
                .fold(f64::INFINITY, f64::min),
        );
    let max_x = data_sim
        .iter()
        .map(|&(x, _)| x)
        .fold(f64::NEG_INFINITY, f64::max)
        .max(
            data_theory
                .iter()
                .map(|&(x, _)| x)
                .fold(f64::NEG_INFINITY, f64::max),
        );
    let min_y = data_sim
        .iter()
        .chain(&data_theory)
        .map(|&(_, y)| y)
        .fold(f64::INFINITY, f64::min);
    let max_y = data_sim
        .iter()
        .chain(&data_theory)
        .map(|&(_, y)| y)
        .fold(f64::NEG_INFINITY, f64::max);

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Intensity Autocorrelation Function g_2(T)",
            ("sans-serif", 64),
        )
        .margin(10)
        .x_label_area_size(80)
        .y_label_area_size(100)
        .build_cartesian_2d(min_x..max_x, min_y..max_y)?;

    chart
        .configure_mesh()
        .x_desc("log10(T) [s]")
        .y_desc("g_2(T)")
        .label_style(("sans-serif", 32))
        .draw()?;

    chart
        .draw_series(LineSeries::new(data_sim, BLUE.stroke_width(4)))?
        .label("Simulated")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    chart
        .draw_series(LineSeries::new(data_theory, GREEN.stroke_width(4)))?
        .label("Theoretical")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart
        .configure_series_labels()
        .label_font(("sans-serif", 32))
        .background_style(&WHITE)
        .draw()?;

    root.present()?;
    Ok(())
}

fn generate_size_plot(path: &Path, result: &DLSResult) -> anyhow::Result<()> {
    let root = BitMapBackend::new(path, (2000, 1000)).into_drawing_area();
    root.fill(&WHITE)?;

    let min_x = result
        .size_nm
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min)
        .min(
            result
                .size_num_nm
                .iter()
                .cloned()
                .fold(f64::INFINITY, f64::min),
        );
    let max_x = result
        .size_nm
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max)
        .max(
            result
                .size_num_nm
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max),
        );
    let min_y = result
        .size_intensity
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min)
        .min(
            result
                .size_num_dist
                .iter()
                .cloned()
                .fold(f64::INFINITY, f64::min),
        );
    let max_y = result
        .size_intensity
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max)
        .max(
            result
                .size_num_dist
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max),
        );

    let mut chart = ChartBuilder::on(&root)
        .caption("Particle Size Distribution", ("sans-serif", 64))
        .margin(10)
        .x_label_area_size(80)
        .y_label_area_size(100)
        .build_cartesian_2d(min_x..max_x, min_y..max_y)?;

    chart
        .configure_mesh()
        .x_desc("Size (nm)")
        .y_desc("Relative Distribution")
        .label_style(("sans-serif", 32))
        .draw()?;

    chart
        .draw_series(LineSeries::new(
            result
                .size_nm
                .iter()
                .zip(&result.size_intensity)
                .map(|(&x, &y)| (x, y)),
            GREEN.stroke_width(4),
        ))?
        .label("Intensity")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart
        .draw_series(LineSeries::new(
            result
                .size_num_nm
                .iter()
                .zip(&result.size_num_dist)
                .map(|(&x, &y)| (x, y)),
            BLUE.stroke_width(4),
        ))?
        .label("Number")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .configure_series_labels()
        .background_style(&WHITE)
        .label_font(("sans-serif", 32))
        .draw()?;

    root.present()?;
    Ok(())
}

pub fn export_csv(
    dir_path: &str,
    _params: &SimulationParams,
    result: &DLSResult,
) -> anyhow::Result<()> {
    use std::path::Path;

    let dir = Path::new(dir_path);

    {
        let default_name = format!(
            "particle_size_distribution_{}.csv",
            Utc::now().format("%Y-%m-%d_%H-%M-%S")
        );
        let path = dir.join(default_name);
        let mut wtr = csv::Writer::from_path(path)?;
        wtr.write_record(&["Size (nm)", "Intensity Distribution", "Number Distribution"])?;
        let max_len = result.size_nm.len().max(result.size_num_nm.len());
        for i in 0..max_len {
            let size_nm = if i < result.size_nm.len() {
                result.size_nm[i].to_string()
            } else {
                "".to_string()
            };
            let intensity = if i < result.size_intensity.len() {
                result.size_intensity[i].to_string()
            } else {
                "".to_string()
            };
            let num_dist = if i < result.size_num_dist.len() {
                result.size_num_dist[i].to_string()
            } else {
                "".to_string()
            };
            wtr.write_record(&[size_nm, intensity, num_dist])?;
        }
        wtr.flush()?;
    }

    {
        let default_name = format!("g2_tau_{}.csv", Utc::now().format("%Y-%m-%d_%H-%M-%S"));
        let path = dir.join(default_name);
        let mut wtr = csv::Writer::from_path(path)?;
        wtr.write_record(&["Tau (s)", "g2 Numeric Noisy", "g2 Theory"])?;
        for i in 0..result.tau.len() {
            wtr.write_record(&[
                result.tau[i].to_string(),
                result.g2_numeric_noisy[i].to_string(),
                result.g2_theory[i].to_string(),
            ])?;
        }
        wtr.flush()?;
    }

    {
        let default_name = format!("g1_tau_{}.csv", Utc::now().format("%Y-%m-%d_%H-%M-%S"));
        let path = dir.join(default_name);
        let mut wtr = csv::Writer::from_path(path)?;
        wtr.write_record(&["Tau (s)", "g1 Numeric Ideal", "g1 Theory"])?;
        for i in 0..result.tau.len() {
            wtr.write_record(&[
                result.tau[i].to_string(),
                result.g1_numeric_ideal[i].to_string(),
                result.g1_theory[i].to_string(),
            ])?;
        }
        wtr.flush()?;
    }

    {
        let default_name = format!("intensity_{}.csv", Utc::now().format("%Y-%m-%d_%H-%M-%S"));
        let path = dir.join(default_name);
        let mut wtr = csv::Writer::from_path(path)?;
        wtr.write_record(&[
            "Time (s)",
            "Intensity Ideal (a.u.)",
            "Intensity Noisy (a.u.)",
        ])?;
        for i in 0..result.time.len() {
            wtr.write_record(&[
                result.time[i].to_string(),
                result.intensity_ideal[i].to_string(),
                result.intensity_noisy[i].to_string(),
            ])?;
        }
        wtr.flush()?;
    }

    {
        let default_name = format!(
            "raw_particle_sizes_{}.csv",
            Utc::now().format("%Y-%m-%d_%H-%M-%S")
        );
        let path = dir.join(default_name);
        let mut wtr = csv::Writer::from_path(path)?;
        wtr.write_record(&["Particle Size (nm)", "Intensity Weight"])?;
        for i in 0..result.raw_sizes_nm.len() {
            let weight = if i < result.intensity_weights.len() {
                result.intensity_weights[i].to_string()
            } else {
                "".to_string()
            };
            wtr.write_record(&[result.raw_sizes_nm[i].to_string(), weight])?;
        }
        wtr.flush()?;
    }

    Ok(())
}
