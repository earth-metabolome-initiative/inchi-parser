//! Tests for invalid InChI strings
use core::num::NonZero;

use inchi_parser::{
    errors::{AtomConnectionTokenError, Error},
    inchi::InChI,
};
use molecular_formulas::errors::ParserError;

#[test]
fn test_missing_inchi_prefix() {
    let inchi_str = "1S/C2H6O/c1-2-3/h3H,2H2,1H3"; // Missing "InChI=" prefix
    let result = inchi_str.parse::<InChI>();
    assert!(matches!(result, Err(Error::MissingInchiPrefix)));
}

#[test]
fn test_missing_version_prefix() {
    let inchi_str = "InChI=/C2H6O/c1-2-3/h3H,2H2,1H3"; // Missing version prefix
    let result = inchi_str.parse::<InChI>();
    assert!(matches!(result, Err(Error::MissingVersionPrefix)));
}

#[test]
fn test_missing_forward_slash() {
    let inchi_str = "InChI=1SC2H6O/c1-2-3/h3H,2H2,1H3"; // Missing '/' after version
    let result = inchi_str.parse::<InChI>();
    assert!(matches!(result, Err(Error::MissingForwardSlash(_))));
}

#[test]
fn test_missing_forward_slash_in_layers() {
    let inchi_str = "InChI=1S/C2H6O"; // Missing '/' before layers
    let result = inchi_str.parse::<InChI>();
    assert!(matches!(result, Err(Error::MissingForwardSlash(_))));
}

#[test]
fn test_not_hill_sorted() {
    let inchi_str = "InChI=1S/C2OH6/"; // Missing '/' before layers
    let result = inchi_str.parse::<InChI>();
    assert!(matches!(result, Err(Error::MolecularFormulaParserError(ParserError::NotHillOrdered))));
}

#[test]
fn test_unbalanced_parenthesis() {
    let inchi_str = "InChI=1S/C2H6O/c1)/";
    let result = inchi_str.parse::<InChI>();
    assert!(matches!(
        result,
        Err(Error::AtomConnectionTokenError(
            AtomConnectionTokenError::ClosingBracketBeforeOpeningBracket
        ))
    ));
}

#[test]
fn test_self_loop_comma() {
    let inchi_str = "InChI=1S/C16H25NS/c1-12(2)13-7-5-8-15(3)9-6-10-16(4,16(13)15)17-11-18/h13-14H,1,5-10H2,2-4H3/t13-,14-,15+,16-/m1/s1";
    let result = inchi_str.parse::<InChI>();
    assert_eq!(
        result,
        Err(Error::AtomConnectionTokenError(AtomConnectionTokenError::SelfLoopDetected(16)))
    );
}

#[test]
fn test_self_loop_dash() {
    let inchi_str = "InChI=1S/C16H25NS/c1-12(2)13-7-5-5-15(3)9-6-10-16(4,16(13)15)17-11-18/h13-14H,1,5-10H2,2-4H3/t13-,14-,15+,16-/m1/s1";
    let result = inchi_str.parse::<InChI>();
    assert_eq!(
        result,
        Err(Error::AtomConnectionTokenError(AtomConnectionTokenError::SelfLoopDetected(5)))
    );
}

#[test]
fn test_self_loop_parenthesis() {
    let inchi_str = "InChI=1S/C16H25NS/c1-12(12)13-7-5-5-15(3)9-6-10-16(4,16(13)15)17-11-18/h13-14H,1,5-10H2,2-4H3/t13-,14-,15+,16-/m1/s1";
    let result = inchi_str.parse::<InChI>();
    assert_eq!(
        result,
        Err(Error::AtomConnectionTokenError(AtomConnectionTokenError::SelfLoopDetected(12)))
    );
}
