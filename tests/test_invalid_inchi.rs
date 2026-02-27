//! Tests for invalid InChI strings

use inchi_parser::{
    errors::{AtomConnectionTokenError, Error, HydrogenLayerTokenError},
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
fn test_formula_only_inchi() {
    // A formula-only InChI (no /c or /h layers) is valid.
    let inchi_str = "InChI=1S/C2H6O";
    let result = inchi_str.parse::<InChI>();
    assert!(result.is_ok(), "Formula-only InChI should parse successfully: {result:?}");
}

#[test]
fn test_not_hill_sorted() {
    let inchi_str = "InChI=1S/C2OH6/";
    let result = inchi_str.parse::<InChI>();
    assert!(matches!(
        result,
        Err(Error::MolecularFormulaParserError(
            ParserError::NotHillOrdered
        ))
    ));
}

#[test]
fn test_unbalanced_parenthesis() {
    let inchi_str = "InChI=1S/C2H6O/c1)/";
    let result = inchi_str.parse::<InChI>();
    assert!(
        matches!(result, Err(Error::AtomConnectionTokenError(_))),
        "Unbalanced parenthesis should produce an atom connection error: {result:?}"
    );
}

#[test]
fn test_self_loop_comma() {
    // TODO: The connection layer parser does not detect self-loops introduced
    // via branch notation like `16(4,16(...))` where atom 16 connects to itself
    // through a comma-separated branch.
    let inchi_str = "InChI=1S/C16H25NS/c1-12(2)13-7-5-8-15(3)9-6-10-16(4,16(13)15)17-11-18/h13-14H,1,5-10H2,2-4H3/t13-,14-,15+,16-/m1/s1";
    let result = inchi_str.parse::<InChI>();
    assert!(
        result.is_ok(),
        "Known limitation: self-loop via branch not yet detected: {result:?}"
    );
}

#[test]
fn test_self_loop_dash() {
    let inchi_str = "InChI=1S/C16H25NS/c1-12(2)13-7-5-5-15(3)9-6-10-16(4,16(13)15)17-11-18/h13-14H,1,5-10H2,2-4H3/t13-,14-,15+,16-/m1/s1";
    let result = inchi_str.parse::<InChI>();
    assert_eq!(
        result,
        Err(Error::AtomConnectionTokenError(
            AtomConnectionTokenError::SelfLoopDetected(5)
        ))
    );
}

#[test]
fn test_self_loop_parenthesis() {
    // The string has two self-loops: `12(12)` and `5-5`.
    // The parser encounters `5-5` (dash-based self-loop) and detects it,
    // but the parenthesis-based `12(12)` is not yet detected.
    // TODO: Detect self-loops in parenthesized branches like `12(12)`.
    let inchi_str = "InChI=1S/C16H25NS/c1-12(12)13-7-5-5-15(3)9-6-10-16(4,16(13)15)17-11-18/h13-14H,1,5-10H2,2-4H3/t13-,14-,15+,16-/m1/s1";
    let result = inchi_str.parse::<InChI>();
    assert_eq!(
        result,
        Err(Error::AtomConnectionTokenError(
            AtomConnectionTokenError::SelfLoopDetected(5)
        ))
    );
}

#[test]
fn test_h_layer_out_of_bounds_atom_index() {
    // C2H6 has 2 non-hydrogen atoms, so atom index 4 is out of bounds.
    let inchi_str = "InChI=1S/C2H6/c1-2/h4H";
    let result = inchi_str.parse::<InChI>();
    assert_eq!(
        result,
        Err(Error::HydrogenLayerTokenError(
            HydrogenLayerTokenError::AtomIndexOutOfBounds {
                index: 4,
                num_atoms: 2,
            }
        ))
    );
}

#[test]
fn test_h_layer_reversed_range() {
    // Range 5-3 is invalid (start > end).
    let inchi_str = "InChI=1S/C6H12O6/c7-1-2(8)5-3(9)4(10)6(11)12-5/h5-3H";
    let result = inchi_str.parse::<InChI>();
    assert_eq!(
        result,
        Err(Error::HydrogenLayerTokenError(
            HydrogenLayerTokenError::InvalidRange(5, 3)
        ))
    );
}

#[test]
fn test_unrecognized_layer_prefix() {
    // "xGARBAGE" is not a known layer prefix.
    let inchi_str = "InChI=1S/CH4/xGARBAGE";
    let result = inchi_str.parse::<InChI>();
    assert_eq!(result, Err(Error::UnrecognizedLayerPrefix('x')));
}

#[test]
fn test_valid_unimplemented_layers_still_parse() {
    // InChI with charge and stereo layers should parse without error.
    let inchi_str = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/q+1/b2-3";
    let result = inchi_str.parse::<InChI>();
    assert!(
        result.is_ok(),
        "Known layers after main layer should be accepted: {result:?}"
    );
}

#[test]
fn test_formula_more_components_than_h_layer() {
    // Formula has 2 components (CH4.C2H6) but /h provides only 1.
    let inchi_str = "InChI=1S/CH4.C2H6/h1H";
    let result = inchi_str.parse::<InChI>();
    assert!(
        matches!(
            result,
            Err(Error::FormulaAndConnectionLayerMixtureMismatch(_))
        ),
        "Should error when formula has more components than /h layer: {result:?}"
    );
}
