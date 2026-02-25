//! Test suite for validating InChI parsing
//! against PubChem data.
//!
//! # Running Tests
//!
//! To run this test (validates all InChI in the PubChem dataset),
//! ensure:
//!
//! ```bash
//! cargo test --release --test test_pubchem_inchi_validation -- --ignored --nocapture
//! ```
use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, Write},
    path::Path,
};

use csv::ReaderBuilder;
use flate2::read::GzDecoder;
use inchi_parser::{errors::Error, inchi::InChI};
use indicatif::{ProgressBar, ProgressStyle};
use serde::Deserialize;

/// Structure representing a PubChem compound entry from the CID-InChI-Key file.
#[derive(Debug, Deserialize)]
struct PubChemCompound {
    /// InChI string
    inchi: String,
}

/// Read and validate PubChem data from the CID-InChI-Key file.
fn validate_pubchem_inchi(file_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let mut csv_reader =
        ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_reader(reader);

    let pb = ProgressBar::new(123_455_803);

    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    let start = std::time::Instant::now();
    let mut parsed_count = 0u64;
    let mut error_examples: HashMap<String, Vec<String>> = HashMap::new();

    for result in csv_reader.deserialize::<PubChemCompound>() {
        let result = result?;
        pb.inc(1);

        let inchi_str = result.inchi.clone();
        match inchi_str.parse::<InChI>() {
            Ok(_) => {}
            Err(Error::UnimplementedFeature(msg)) => {}
            Err(e) => {
                let error_key = e.to_string();

                let entry = error_examples.entry(error_key).or_default();

                if entry.len() < 2 {
                    entry.push(inchi_str);
                }
            }
        }

        parsed_count += 1;
    }

    let time_required = start.elapsed().as_secs_f64();
    #[allow(clippy::cast_precision_loss)]
    let time_per_compound = time_required / parsed_count as f64;

    pb.finish_with_message("Validation complete");

    println!(
        "Time taken: {:.2} seconds ({:.6} milliseconds per compound)",
        time_required,
        time_per_compound * 1000.0
    );
    println!("Parsed compounds: {parsed_count}");

    if !error_examples.is_empty() {
        let mut file = File::create("tests/failed_inchis.txt")?;

        for (error, examples) in &error_examples {
            writeln!(file, "Error: {error}")?;
            for example in examples {
                writeln!(file, "  Example: {example}")?;
            }
            writeln!(file)?;
        }

        println!("Wrote {} error types to failed_inchis.txt", error_examples.len());
    }

    Ok(())
}

#[test]
#[ignore = "This test requires the sorted_pubchem_inchi.tsv.gz file to be present and is time-consuming."]
/// Validate InChI parsing against PubChem
/// CID-InChI-Key data.
///
/// The document, weighing compressed approximately 6.79 GB, can be downloaded
/// from:
///
/// <https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz>
fn test_pubchem_validation() {
    let file_path = Path::new("tests/sorted_pubchem_inchi.tsv.gz");

    if !file_path.exists() {
        eprintln!("sorted_pubchem_inchi.tsv.gz file not found. Skipping test.");
        return;
    }

    println!("Validating all PubChem InChI formulas...");

    validate_pubchem_inchi(file_path).expect("Validation failed");
}
