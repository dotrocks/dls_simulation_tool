use crate::structs::{DLSResult, SimulationParams};
use genpdf::Element;
use genpdf::{Alignment, elements, style};
use std::fs::File;

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

    let file = File::create(path)?;
    doc.render(file)?;

    Ok(())
}
