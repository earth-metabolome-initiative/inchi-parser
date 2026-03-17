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
    io::{BufRead, BufReader, Write},
    path::Path,
};

use flate2::read::GzDecoder;
use inchi_parser::inchi::InChI;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

/// The URL for the PubChem CID-InChI-Key gzipped file.
const PUBCHEM_URL: &str = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz";

/// Local file name for the downloaded PubChem data.
const PUBCHEM_FILE: &str = "tests/CID-InChI-Key.gz";

/// Local file name for the decompressed PubChem data.
const PUBCHEM_TSV_FILE: &str = "tests/CID-InChI-Key.tsv";

/// Estimated number of compounds in the PubChem dataset.
const ESTIMATED_COMPOUNDS: u64 = 123_000_000;

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

/// Decompress the gzipped PubChem file to a plain TSV if it doesn't already
/// exist.
fn ensure_decompressed(gz_path: &Path, tsv_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    if tsv_path.exists() {
        println!("Decompressed file already exists at {}", tsv_path.display());
        return Ok(());
    }

    println!("Decompressing {} to {}...", gz_path.display(), tsv_path.display());

    let gz_file = File::open(gz_path)?;
    let decoder = GzDecoder::new(gz_file);
    let mut reader = BufReader::new(decoder);
    let mut writer = std::io::BufWriter::new(File::create(tsv_path)?);

    let mut buf = vec![0u8; 64 * 1024];
    let mut total_bytes: u64 = 0;
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {bytes} decompressed")
            .unwrap(),
    );

    loop {
        let n = std::io::Read::read(&mut reader, &mut buf)?;
        if n == 0 {
            break;
        }
        std::io::Write::write_all(&mut writer, &buf[..n])?;
        total_bytes += n as u64;
        pb.set_position(total_bytes);
    }

    pb.finish_with_message("Decompression complete");
    println!("Decompressed {total_bytes} bytes to {}", tsv_path.display());

    Ok(())
}

/// Download (if needed) and decompress (if needed) the PubChem TSV file.
fn ensure_pubchem_tsv() {
    let gz_path = Path::new(PUBCHEM_FILE);
    let tsv_path = Path::new(PUBCHEM_TSV_FILE);

    let rt = tokio::runtime::Runtime::new().expect("Failed to create tokio runtime");
    rt.block_on(ensure_pubchem_file(gz_path)).expect("Failed to ensure PubChem file is available");

    ensure_decompressed(gz_path, tsv_path).expect("Failed to decompress PubChem file");
}

/// Extract the InChI string (second tab-separated field) from a TSV line.
fn extract_inchi(line: &str) -> Option<&str> {
    line.split('\t').nth(1)
}

/// Create a progress bar with the standard style for PubChem iteration.
fn make_progress_bar() -> ProgressBar {
    let pb = ProgressBar::new(ESTIMATED_COMPOUNDS);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );
    pb
}

/// Number of lines to read per chunk for parallel processing.
const CHUNK_SIZE: usize = 100_000;

/// Read the next chunk of lines from a buffered reader into the provided
/// buffer.
///
/// Returns the number of lines read (0 at EOF).
fn read_chunk(reader: &mut impl BufRead, chunk: &mut Vec<String>) -> std::io::Result<usize> {
    chunk.clear();
    let mut line = String::new();
    for _ in 0..CHUNK_SIZE {
        line.clear();
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        chunk.push(std::mem::take(&mut line));
    }
    Ok(chunk.len())
}

/// Read and validate PubChem data from the decompressed TSV file.
fn validate_pubchem_inchi(file_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let mut reader = BufReader::new(file);

    let pb = make_progress_bar();
    let mut parsed_count = 0u64;
    let mut error_examples: HashMap<String, Vec<String>> = HashMap::new();

    let start = std::time::Instant::now();
    let mut chunk = Vec::with_capacity(CHUNK_SIZE);

    loop {
        let n = read_chunk(&mut reader, &mut chunk)?;
        if n == 0 {
            break;
        }

        // Process chunk in parallel, collect errors
        let chunk_errors: Vec<(String, String)> = chunk
            .par_iter()
            .filter_map(|line| {
                let inchi = extract_inchi(line)?;
                if let Err(e) = inchi.parse::<InChI>() {
                    Some((e.to_string(), inchi.to_string()))
                } else {
                    None
                }
            })
            .collect();

        // Merge errors sequentially (cheap — errors are rare)
        for (error_key, inchi) in chunk_errors {
            let entry = error_examples.entry(error_key).or_default();
            if entry.len() < 2 {
                entry.push(inchi);
            }
        }

        parsed_count += n as u64;
        pb.set_position(parsed_count);
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
#[ignore = "This test downloads a ~6.79 GB file and is time-consuming."]
/// Validate InChI parsing against PubChem CID-InChI-Key data.
///
/// The document is automatically downloaded from:
///
/// <https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz>
fn test_pubchem_validation() {
    ensure_pubchem_tsv();
    println!("Validating all PubChem InChI formulas...");
    validate_pubchem_inchi(Path::new(PUBCHEM_TSV_FILE)).expect("Validation failed");
}

/// Bucket statistics for a single category of display difference.
struct BucketStats {
    /// Number of InChIs with this difference.
    count: u64,
    /// Up to 3 example pairs of (original_segment, displayed_segment).
    examples: Vec<(String, String)>,
}

/// Layer prefix to human-readable bucket name.
fn layer_name(prefix: &str) -> &'static str {
    match prefix {
        "formula" => "formula",
        "c" => "connection_layer",
        "h" => "hydrogen_layer",
        "q" => "charge_layer",
        "p" => "proton_layer",
        "b" => "double_bond_stereo",
        "t" => "tetrahedral_stereo",
        "m" => "allene_stereo",
        "s" => "stereo_info",
        "i" => "isotope_layer",
        _ => "unknown_layer",
    }
}

/// Known layer prefix characters in InChI order.
const LAYER_PREFIXES: &[&str] = &["c", "h", "q", "p", "b", "t", "m", "s", "i"];

/// Split an InChI string (after the `InChI=1S/` prefix) into a map from
/// layer prefix to segment content, plus the formula as key `"formula"`.
fn split_into_layers(body: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    let segments: Vec<&str> = body.split('/').collect();
    if segments.is_empty() {
        return map;
    }

    // First segment is the formula (unless it starts with a known layer
    // prefix, which would mean no formula — e.g. proton-only).
    let first = segments[0];
    if !first.is_empty()
        && !LAYER_PREFIXES.iter().any(|p| first.starts_with(*p) && first.len() > p.len())
    {
        // Might be a formula or might be proton-only starting with 'p'
        // Heuristic: formulas start with an uppercase letter
        if first.starts_with(|c: char| c.is_ascii_uppercase()) {
            map.insert("formula".to_string(), first.to_string());
        }
    }

    // If formula wasn't inserted, the first segment might be a layer
    let start_idx = usize::from(map.contains_key("formula"));

    for &seg in &segments[start_idx..] {
        if seg.is_empty() {
            continue;
        }
        // Find the prefix: it's the leading alphabetic character(s) that match
        // a known prefix, but for isotope /i the sublayer /h follows inside
        // the isotope layer. We just use the first char.
        let prefix = &seg[..1];
        // For the isotope hydrogen sublayer `hD2T3` that follows `/i`, we
        // need special handling, but since we split on `/` each segment is
        // independent. The isotope `/h` sublayer will show up as a segment
        // starting with `h` — which conflicts with the main hydrogen layer.
        // InChI ordering guarantees main `/h` comes before `/i`, so the
        // second `h` segment is the isotope hydrogen sublayer. We handle
        // this by appending to existing key with a suffix.
        let key = prefix.to_string();
        if map.contains_key(&key) {
            if key == "h" && map.contains_key("i") {
                // Second /h after /i is the isotope hydrogen sublayer
                map.insert("h_isotope".to_string(), seg.to_string());
            }
            // Otherwise ignore — duplicate arose from comparison of
            // different n* compression forms, not a real second layer.
        } else {
            map.insert(key, seg.to_string());
        }
    }
    map
}

/// Result of comparing a single InChI's Display output against its original
/// string.
enum DisplayResult {
    /// InChI failed to parse.
    ParseError,
    /// Display output matches the original exactly.
    ExactMatch,
    /// Display output differs; contains the differing layer buckets.
    Difference(Vec<(String, String, String)>),
}

/// Compare display output against PubChem originals and produce a bucket
/// report.
fn compare_pubchem_display(file_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let mut reader = BufReader::new(file);

    let pb = make_progress_bar();
    let mut total_parsed: u64 = 0;
    let mut exact_matches: u64 = 0;
    let mut parse_errors: u64 = 0;
    let mut with_differences: u64 = 0;
    let mut buckets: HashMap<String, BucketStats> = HashMap::new();

    let start = std::time::Instant::now();
    let mut chunk = Vec::with_capacity(CHUNK_SIZE);
    let mut lines_read = 0u64;

    loop {
        let n = read_chunk(&mut reader, &mut chunk)?;
        if n == 0 {
            break;
        }

        // Process chunk in parallel
        let results: Vec<DisplayResult> = chunk
            .par_iter()
            .filter_map(|line| {
                let inchi_str = extract_inchi(line)?;

                let Ok(parsed) = inchi_str.parse::<InChI>() else {
                    return Some(DisplayResult::ParseError);
                };

                let displayed = parsed.to_string();
                if displayed == inchi_str {
                    return Some(DisplayResult::ExactMatch);
                }

                // Collect per-layer diffs
                let orig_body = inchi_str
                    .strip_prefix("InChI=1S/")
                    .or_else(|| inchi_str.strip_prefix("InChI=1/"))
                    .unwrap_or(inchi_str);
                let disp_body = displayed
                    .strip_prefix("InChI=1S/")
                    .or_else(|| displayed.strip_prefix("InChI=1/"))
                    .unwrap_or(&displayed);

                let orig_layers = split_into_layers(orig_body);
                let disp_layers = split_into_layers(disp_body);

                let mut all_keys: Vec<String> = orig_layers.keys().cloned().collect();
                for k in disp_layers.keys() {
                    if !all_keys.contains(k) {
                        all_keys.push(k.clone());
                    }
                }

                let mut diffs = Vec::new();
                for key in &all_keys {
                    let orig_val = orig_layers.get(key);
                    let disp_val = disp_layers.get(key);
                    if orig_val != disp_val {
                        let bucket_name =
                            if key == "h_isotope" { "isotope_h_sublayer" } else { layer_name(key) };
                        let orig_str =
                            orig_val.map_or_else(|| "(missing)".to_string(), Clone::clone);
                        let disp_str =
                            disp_val.map_or_else(|| "(missing)".to_string(), Clone::clone);
                        diffs.push((bucket_name.to_string(), orig_str, disp_str));
                    }
                }

                // Capture uncategorized diffs where no individual layer differed
                if diffs.is_empty() {
                    diffs.push(("uncategorized".to_string(), inchi_str.to_string(), displayed));
                }

                Some(DisplayResult::Difference(diffs))
            })
            .collect();

        // Merge results sequentially
        for result in results {
            match result {
                DisplayResult::ParseError => parse_errors += 1,
                DisplayResult::ExactMatch => {
                    total_parsed += 1;
                    exact_matches += 1;
                }
                DisplayResult::Difference(diffs) => {
                    total_parsed += 1;
                    with_differences += 1;
                    for (bucket_name, orig_str, disp_str) in diffs {
                        let max_examples = if bucket_name == "uncategorized" { 20 } else { 3 };
                        let stats = buckets
                            .entry(bucket_name)
                            .or_insert_with(|| BucketStats { count: 0, examples: Vec::new() });
                        stats.count += 1;
                        if stats.examples.len() < max_examples {
                            stats.examples.push((orig_str, disp_str));
                        }
                    }
                }
            }
        }

        lines_read += n as u64;
        pb.set_position(lines_read);
    }

    let time_required = start.elapsed().as_secs_f64();

    pb.finish_with_message("Comparison complete");

    println!("Time taken: {time_required:.2} seconds");

    write_comparison_report(total_parsed, exact_matches, with_differences, parse_errors, &buckets)
}

/// Write the comparison report to stdout and to `tests/display_comparison.txt`.
fn write_comparison_report(
    total_parsed: u64,
    exact_matches: u64,
    with_differences: u64,
    parse_errors: u64,
    buckets: &HashMap<String, BucketStats>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut sorted_buckets: Vec<(&String, &BucketStats)> = buckets.iter().collect();
    sorted_buckets.sort_by_key(|entry| std::cmp::Reverse(entry.1.count));

    #[allow(clippy::cast_precision_loss)]
    let exact_pct =
        if total_parsed > 0 { exact_matches as f64 / total_parsed as f64 * 100.0 } else { 0.0 };
    #[allow(clippy::cast_precision_loss)]
    let diff_pct =
        if total_parsed > 0 { with_differences as f64 / total_parsed as f64 * 100.0 } else { 0.0 };

    println!("\n=== PubChem Display Comparison ===");
    println!("Total parsed:        {total_parsed}");
    println!("Exact matches:       {exact_matches} ({exact_pct:.1}%)");
    println!("With differences:    {with_differences} ({diff_pct:.1}%)");
    println!("Parse errors:        {parse_errors}");
    println!();
    println!("Difference buckets (sorted by frequency):");

    for (name, stats) in &sorted_buckets {
        #[allow(clippy::cast_precision_loss)]
        let pct =
            if total_parsed > 0 { stats.count as f64 / total_parsed as f64 * 100.0 } else { 0.0 };
        println!("  {name:<25} {:<10} ({pct:.1}%)", stats.count);
    }

    let mut out = File::create("tests/display_comparison.txt")?;
    writeln!(out, "=== PubChem Display Comparison — Detailed Examples ===")?;
    writeln!(out, "Total parsed:      {total_parsed}")?;
    writeln!(out, "Exact matches:     {exact_matches} ({exact_pct:.1}%)")?;
    writeln!(out, "With differences:  {with_differences} ({diff_pct:.1}%)")?;
    writeln!(out, "Parse errors:      {parse_errors}")?;
    writeln!(out)?;

    for (name, stats) in &sorted_buckets {
        #[allow(clippy::cast_precision_loss)]
        let pct =
            if total_parsed > 0 { stats.count as f64 / total_parsed as f64 * 100.0 } else { 0.0 };
        writeln!(out, "--- {name}: {} ({pct:.1}%) ---", stats.count)?;
        for (orig, disp) in &stats.examples {
            writeln!(out, "  original:  {orig}")?;
            writeln!(out, "  displayed: {disp}")?;
            writeln!(out)?;
        }
    }

    println!("\nExamples per bucket written to tests/display_comparison.txt");

    Ok(())
}

#[test]
#[ignore = "This test downloads a ~6.79 GB file and is time-consuming."]
/// Compare Display output to original PubChem InChI strings and classify
/// differences.
fn test_pubchem_display_comparison() {
    ensure_pubchem_tsv();
    println!("Comparing Display output against PubChem InChI strings...");
    compare_pubchem_display(Path::new(PUBCHEM_TSV_FILE)).expect("Comparison failed");
}
