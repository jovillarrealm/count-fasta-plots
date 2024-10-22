use plotters::prelude::*;
use plotters::coord::types::RangedCoordf64;
use serde::Deserialize;
use std::error::Error;
use statistical::median;

#[derive(Debug, Deserialize)]
struct GenomeStats {
    assembly_length: f64,
    number_of_sequences: f64,
    #[serde(rename = "N50")]
    n50: f64,
    #[serde(rename = "GC_percentage")]
    gc_percentage: f64,
    #[serde(rename = "N_percentage")]
    n_percentage: f64,
}

fn create_boxplot<DB: DrawingBackend>(
    plot: &mut ChartContext<DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    data: &[f64],
    y_position: f64,
) -> Result<(), Box<dyn Error>> where <DB as plotters::prelude::DrawingBackend>::ErrorType: 'static {
    // Calculate box plot statistics
    let mut sorted_data = data.to_vec();
    sorted_data.sort_by(|a, b| a.partial_cmp(b).unwrap());
    
    let q1_idx = (sorted_data.len() as f64 * 0.25) as usize;
    let q3_idx = (sorted_data.len() as f64 * 0.75) as usize;
    
    let q1 = sorted_data[q1_idx];
    let q3 = sorted_data[q3_idx];
    let med = median(&sorted_data);
    
    // Calculate IQR and whisker bounds
    let iqr = q3 - q1;
    let lower_bound = q1 - 1.5 * iqr;
    let upper_bound = q3 + 1.5 * iqr;
    
    // Find actual whisker ends (last non-outlier points)
    let whisker_min = sorted_data.iter()
        .find(|&&x| x >= lower_bound)
        .copied()
        .unwrap_or(q1);
    let whisker_max = sorted_data.iter()
        .rev()
        .find(|&&x| x <= upper_bound)
        .copied()
        .unwrap_or(q3);
    
    // Draw box
    plot.draw_series(std::iter::once(Rectangle::new(
        [(q1, y_position - 0.3), (q3, y_position + 0.3)],
        BLUE.mix(0.3),
    )))?;
    
    // Draw median line
    plot.draw_series(std::iter::once(Rectangle::new(
        [(med - 0.1, y_position - 0.3), (med + 0.1, y_position + 0.3)],
        RED,
    )))?;
    
    // Draw whiskers
    plot.draw_series(std::iter::once(PathElement::new(
        vec![(whisker_min, y_position), (q1, y_position)],
        BLACK,
    )))?;
    plot.draw_series(std::iter::once(PathElement::new(
        vec![(q3, y_position), (whisker_max, y_position)],
        BLACK,
    )))?;
    
    // Draw outlier points
    let outliers: Vec<_> = data.iter()
        .filter(|&&x| x < lower_bound || x > upper_bound)
        .collect();
    
    plot.draw_series(outliers.iter().map(|&&x| {
        Circle::new((x, y_position), 3, BLACK.filled())
    }))?;
    
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut data = Vec::new();
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b';')
        .from_path(std::env::args().nth(1).expect("Please provide a CSV file"))?;
    
    for result in rdr.deserialize() {
        let record: GenomeStats = result?;
        data.push(record);
    }

    let root = BitMapBackend::new("count-fasta.png", (800, 1120))
        .into_drawing_area();
    root.fill(&WHITE)?;

    let plots = root.split_evenly((5, 1));
    
    let metrics: Vec<(&str, &str, Box<dyn Fn(&GenomeStats) -> f64>)> = vec![
        ("Assembly size (bp.)", "bp.", Box::new(|x| x.assembly_length)),
        ("Scaffold count", "Count", Box::new(|x| x.number_of_sequences)),
        ("N50 (bp.)", "bp.", Box::new(|x| x.n50)),
        ("GC ratio (%)", "GC ratio (%)", Box::new(|x| x.gc_percentage)),
        ("N's ratio (%)", "Ratio (%)", Box::new(|x| x.n_percentage)),
    ];

    for (i, (plot_area, (title, xlabel, metric))) in plots.iter().zip(metrics.iter()).enumerate() {
        let values: Vec<f64> = data.iter().map(metric).collect();
        let min = values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max = values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let padding = (max - min) * 0.1;

        let mut chart = ChartBuilder::on(plot_area)
            .margin(5)
            .caption(title, ("sans-serif", 20))
            .set_label_area_size(LabelAreaPosition::Left, 40)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .build_cartesian_2d(min - padding..max + padding, 0.0..2.0)?;

        chart
            .configure_mesh()
            .disable_y_mesh()
            .y_desc(format!("{}", (b'A' + i as u8) as char))
            .x_desc(xlabel.to_string())
            .draw()?;

        create_boxplot(&mut chart, &values, 1.0)?;
    }

    Ok(())
}