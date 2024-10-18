use std::error::Error;
use csv::ReaderBuilder;
use plotters::prelude::*;
use serde::Deserialize;
use plotters::prelude::full_palette::BLACK;
use plotters::coord::Shift;

#[derive(Debug, Deserialize)]
struct FastaStats {
    assembly_length: f64,
    number_of_sequences: f64,
    #[serde(rename = "N50")]
    n50: f64,
    #[serde(rename = "GC_percentage")]
    gc_percentage: f64,
    #[serde(rename = "N_percentage")]
    n_percentage: f64,
}

fn create_box_plot<DB: DrawingBackend>(
    ctx: &DrawingArea<DB, Shift>,
    data: &[f64],
    title: &str,
    x_label: &str,
    y_label: &str,
) -> Result<(), Box<dyn Error>>
where
    DB::ErrorType: 'static,
{
    let min = data.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = max - min;
    let padding = range * 0.1;

    let mut chart = ChartBuilder::on(ctx)
        .margin(5)
        .caption(title, ("sans-serif", 20))
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d((min - padding)..(max + padding), 0f64..1f64)?;

    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .draw()?;

    // Calculate box plot statistics
    let mut sorted_data: Vec<f64> = data.to_vec();
    sorted_data.sort_by(|a, b| a.partial_cmp(b).unwrap());
    
    let q1_idx = (sorted_data.len() as f64 * 0.25) as usize;
    let q2_idx = (sorted_data.len() as f64 * 0.5) as usize;
    let q3_idx = (sorted_data.len() as f64 * 0.75) as usize;
    
    let q1 = sorted_data[q1_idx];
    let median = sorted_data[q2_idx];
    let q3 = sorted_data[q3_idx];
    
    let iqr = q3 - q1;
    let lower_fence = q1 - 1.5 * iqr;
    let upper_fence = q3 + 1.5 * iqr;

    // Draw box
    chart.draw_series(std::iter::once(Rectangle::new(
        [(q1, 0.25), (q3, 0.75)],
        BLUE.mix(0.3).filled(),
    )))?;

    // Draw median line
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(median, 0.25), (median, 0.75)],
        RED.stroke_width(2),
    )))?;

    // Draw whiskers
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(lower_fence, 0.5), (q1, 0.5)],
        BLACK.stroke_width(1),
    )))?;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(q3, 0.5), (upper_fence, 0.5)],
        BLACK.stroke_width(1),
    )))?;

    // Draw outliers
    let outliers: Vec<(f64, f64)> = sorted_data
        .iter()
        .filter(|&&x| x < lower_fence || x > upper_fence)
        .map(|&x| (x, 0.5))
        .collect();

    chart.draw_series(outliers.iter().map(|point| {
        Circle::new(*point, 3, BLACK.mix(0.5).filled())
    }))?;

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let filename = std::env::args().nth(1).expect("Please provide a CSV file path.\nOutput in the shape of the output of count-fasta-rs");
    
    // Read CSV file
    let mut reader = ReaderBuilder::new()
        .delimiter(b';')
        .from_path(&filename)?;
    
    let records: Vec<FastaStats> = reader
        .deserialize()
        .collect::<Result<Vec<FastaStats>, _>>()?;

    // Prepare data vectors
    let assembly_lengths: Vec<f64> = records.iter().map(|r| r.assembly_length).collect();
    let sequence_counts: Vec<f64> = records.iter().map(|r| r.number_of_sequences).collect();
    let n50s: Vec<f64> = records.iter().map(|r| r.n50).collect();
    let gc_percentages: Vec<f64> = records.iter().map(|r| r.gc_percentage).collect();
    let n_percentages: Vec<f64> = records.iter().map(|r| r.n_percentage).collect();

    // Create the output file
    let root = BitMapBackend::new("count-fasta.png", (800, 1120))
        .into_drawing_area();
    root.fill(&WHITE)?;

    // Split the drawing area into 5 parts
    let areas = root.split_evenly((5, 1));

    // Create each box plot
    create_box_plot(
        &areas[0],
        &assembly_lengths,
        "Assembly size (bp.)",
        "bp.",
        "A",
    )?;
    create_box_plot(
        &areas[1],
        &sequence_counts,
        "Scaffold count",
        "Count",
        "B",
    )?;
    create_box_plot(
        &areas[2],
        &n50s,
        "N50 (bp.)",
        "bp.",
        "C",
    )?;
    create_box_plot(
        &areas[3],
        &gc_percentages,
        "GC ratio (%)",
        "GC ratio (%)",
        "D",
    )?;
    create_box_plot(
        &areas[4],
        &n_percentages,
        "N's ratio (%)",
        "Ratio (%)",
        "E",
    )?;

    root.present()?;
    println!("Plot has been saved as 'count-fasta.png'");

    Ok(())
}