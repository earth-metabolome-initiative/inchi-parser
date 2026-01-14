//! Tests for invalid InChI strings
use inchi_parser::errors::Errors;
use inchi_parser::inchi::InChI;

#[test]
fn test_missing_inchi_prefix() {
    let inchi_str = "1S/C2H6O/c1-2-3/h3H,2H2,1H3"; // Missing "InChI=" prefix
    let result = inchi_str.parse::<InChI>();
    assert!(matches!(result, Err(Errors::MissingInchiPrefix)));
}

#[test]
fn test_missing_version_prefix() {
    let inchi_str = "InChI=/C2H6O/c1-2-3/h3H,2H2,1H3"; // Missing version prefix
    let result = inchi_str.parse::<InChI>();
    assert!(matches!(result, Err(Errors::MissingVersionPrefix)));
}

#[test]
fn test_missing_forward_slash() {
    let inchi_str = "InChI=1SC2H6O/c1-2-3/h3H,2H2,1H3"; // Missing '/' after version
    let result = inchi_str.parse::<InChI>();
    assert!(matches!(result, Err(Errors::MissingForwardSlash)));
}

#[test]
fn test_missing_forward_slash_in_layers() {
    let inchi_str = "InChI=1S/C2H6O"; // Missing '/' before layers
    let result = inchi_str.parse::<InChI>();
    assert!(matches!(result, Err(Errors::MissingForwardSlash)));
}

#[test]
fn test_not_hill_sorted() {
    let inchi_str = "InChI=1S/C2OH6/"; // Missing '/' before layers
    let result = inchi_str.parse::<InChI>();
    assert!(matches!(result, Err(Errors::NotHillSorted)));
}
