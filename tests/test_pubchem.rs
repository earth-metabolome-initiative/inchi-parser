//! Test suite for validating InChI parsing
//! against PubChem data.
//!
//! # Running Tests
//!
//! To run this test (validates all InChI in the PubChem dataset),
//! ensure:
//!
//! ```bash
//! cargo test --release --test test_pubchem -- --ignored --nocapture
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

/// The URL for the PubChem CID-InChI-Key gzipped file.
const PUBCHEM_URL: &str =
    "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz";

/// Local file name for the downloaded PubChem data.
const PUBCHEM_FILE: &str = "tests/CID-InChI-Key.gz";

/// Structure representing a PubChem compound entry from the CID-InChI-Key file.
///
/// The file is a tab-separated values file with three columns:
/// CID, InChI, and InChIKey.
#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct PubChemCompound {
    /// PubChem Compound ID
    cid: u64,
    /// InChI string
    inchi: String,
    /// InChI Key
    inchi_key: String,
}

/// Download the PubChem CID-InChI-Key file if it doesn't already exist.
async fn ensure_pubchem_file(file_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    if file_path.exists() {
        println!("PubChem file already exists at {}", file_path.display());
        return Ok(());
    }

    println!("Downloading PubChem CID-InChI-Key file from {PUBCHEM_URL}...");
    println!("This file is approximately 6.79 GB and may take a while to download.");

    let response = reqwest::get(PUBCHEM_URL).await?;

    if !response.status().is_success() {
        return Err(format!("HTTP request failed with status: {}", response.status()).into());
    }

    let total_size = response.content_length().unwrap_or(0);

    let pb = ProgressBar::new(total_size);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({bytes_per_sec}, {eta})")
            .unwrap()
            .progress_chars("#>-"),
    );

    let bytes = response.bytes().await?;
    pb.set_position(bytes.len() as u64);

    tokio::fs::write(file_path, &bytes).await?;

    pb.finish_with_message("Download complete");
    println!("Downloaded {} bytes.", bytes.len());

    Ok(())
}

/// Read and validate PubChem data from the CID-InChI-Key file.
fn validate_pubchem_inchi(file_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let mut csv_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(reader);

    // Estimated count based on recent PubChem data
    let pb = ProgressBar::new(123_000_000);

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
        let compound = result?;
        pb.inc(1);

        match compound.inchi.parse::<InChI>() {
            Ok(_) | Err(Error::UnimplementedFeature(_)) => {}
            Err(e) => {
                let error_key = e.to_string();
                let entry = error_examples.entry(error_key).or_default();
                if entry.len() < 2 {
                    entry.push(compound.inchi);
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

        println!(
            "Wrote {} error types to failed_inchis.txt",
            error_examples.len()
        );
    }

    Ok(())
}

#[test]
#[ignore = "This test downloads a ~6.79 GB file and is time-consuming."]
/// Validate InChI parsing against PubChem CID-InChI-Key data.
///
/// The document is automatically downloaded from:
///
/// <https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz>
fn test_pubchem_validation() {
    let file_path = Path::new(PUBCHEM_FILE);

    let rt = tokio::runtime::Runtime::new().expect("Failed to create tokio runtime");
    rt.block_on(ensure_pubchem_file(file_path))
        .expect("Failed to ensure PubChem file is available");

    println!("Validating all PubChem InChI formulas...");

    validate_pubchem_inchi(file_path).expect("Validation failed");
}
