mod gui;
mod simulation;
mod structs;

use crate::gui::DLSApp;

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1500.0, 1000.0])
            .with_min_inner_size([1200.0, 800.0]),
        ..Default::default()
    };

    eframe::run_native(
        "Dynamic Light Scattering Simulator",
        options,
        Box::new(|cc| Ok(Box::new(DLSApp::new(cc)))),
    )
}
