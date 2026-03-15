use alloc::{vec, vec::Vec};

use elements_rs::isotopes::HydrogenIsotope;
use molecular_formulas::{BaselineDigit, InChIFormula, MolecularFormula, try_fold_number};

use crate::{
    errors::Error,
    impls::charge_layer::parse_charge,
    inchi::isotope_layer::{IsotopeAtom, IsotopeComponent, IsotopeHydrogen, IsotopeLayer},
    traits::{
        parse::{FromStrWithContext, PrefixFromStrWithContext},
        prefix::Prefix,
    },
};

/// Parses hydrogen isotope specs from a segment like `hD2`, `hT`, `hDT3`.
fn parse_h_isotope_segment(s: &str) -> Result<Vec<IsotopeHydrogen>, Error<u16>> {
    let s = s.strip_prefix('h').ok_or(Error::InvalidIsotopeValue('h'))?;
    let mut hydrogens = Vec::new();
    let mut chars = s.chars().peekable();

    while chars.peek().is_some() {
        let isotope = match chars.next() {
            Some('D') => HydrogenIsotope::D,
            Some('T') => HydrogenIsotope::T,
            Some('H') => HydrogenIsotope::H1,
            Some(c) => return Err(Error::InvalidIsotopeValue(c)),
            None => unreachable!(),
        };
        let count = try_fold_number::<u16, BaselineDigit, _>(&mut chars)
            .transpose()
            .map_err(|_| Error::InvalidIsotopeValue('0'))?
            .unwrap_or(1);
        hydrogens.push(IsotopeHydrogen { isotope, count });
    }

    Ok(hydrogens)
}

/// Parses atom isotope specs from a string like `1+1,3+2` or `1-1`.
fn parse_atom_specs(s: &str) -> Result<Vec<IsotopeAtom>, Error<u16>> {
    if s.is_empty() {
        return Ok(Vec::new());
    }

    let mut atoms = Vec::new();
    for spec in s.split(',') {
        // Find the position of the sign character (+/-) that separates index from
        // shift.
        let sign_pos = spec.find(['+', '-']).ok_or_else(|| {
            spec.chars().next().map_or(Error::InvalidIsotopeValue('?'), Error::InvalidIsotopeValue)
        })?;

        let index_str = &spec[..sign_pos];
        let shift_str = &spec[sign_pos..];

        // Parse 1-based atom index
        let mut chars = index_str.chars().peekable();
        let one_based = match try_fold_number::<u16, BaselineDigit, _>(&mut chars) {
            Some(Ok(n)) if n > 0 => n,
            _ => {
                return Err(Error::InvalidIsotopeValue(index_str.chars().next().unwrap_or('?')));
            }
        };

        if let Some(&c) = chars.peek() {
            return Err(Error::InvalidIsotopeValue(c));
        }

        let atom_index = one_based - 1;

        let mass_shift = parse_charge(shift_str)?;

        atoms.push(IsotopeAtom { atom_index, mass_shift });
    }

    Ok(atoms)
}

impl FromStrWithContext for IsotopeLayer {
    type Context<'a> = (&'a InChIFormula, Option<&'a str>);
    type Input<'a> = &'a str;
    type Idx = u16;

    fn from_str_with_context(
        input: Self::Input<'_>,
        context: Self::Context<'_>,
    ) -> Result<Self, Error<Self::Idx>> {
        let s = input.strip_prefix(Self::PREFIX).ok_or(Error::WrongPrefix)?;
        let (formula, h_segment) = context;

        // Parse hydrogen isotope specs from the /h sublayer (shared across components
        // when it appears as a single segment without semicolons).
        let shared_hydrogens =
            if let Some(h_seg) = h_segment { parse_h_isotope_segment(h_seg)? } else { Vec::new() };

        // If the atom spec body is empty, build a single-component layer with just H
        // isotopes.
        if s.is_empty() {
            let components =
                vec![IsotopeComponent { atoms: Vec::new(), hydrogens: shared_hydrogens }];
            return Ok(IsotopeLayer { components });
        }

        let mut subformulas = formula.subformulas();
        let mut components = Vec::with_capacity(formula.number_of_mixtures());

        for component_str in s.split(';') {
            // Detect optional n* repetition prefix.
            let digit_count = component_str.bytes().take_while(u8::is_ascii_digit).count();
            let (reps, component_str) =
                if digit_count > 0 && component_str.as_bytes().get(digit_count) == Some(&b'*') {
                    component_str[..digit_count]
                        .parse::<u32>()
                        .ok()
                        .map_or((1u32, component_str), |n| (n, &component_str[digit_count + 1..]))
                } else {
                    (1u32, component_str)
                };

            subformulas.next().ok_or(Error::FormulaAndConnectionLayerMixtureMismatch(
                formula.number_of_mixtures(),
            ))?;

            let atoms = parse_atom_specs(component_str)?;
            let component = IsotopeComponent { atoms, hydrogens: shared_hydrogens.clone() };

            components.push(component.clone());
            for _ in 1..reps {
                subformulas.next().ok_or(Error::FormulaAndConnectionLayerMixtureMismatch(
                    formula.number_of_mixtures(),
                ))?;
                components.push(component.clone());
            }
        }

        if subformulas.next().is_some() {
            return Err(Error::FormulaAndConnectionLayerMixtureMismatch(
                formula.number_of_mixtures(),
            ));
        }

        Ok(IsotopeLayer { components })
    }
}

impl PrefixFromStrWithContext for IsotopeLayer {
    fn try_build_layer(
        input: &mut &str,
        context: Self::Context<'_>,
    ) -> Result<Option<Self>, Error<Self::Idx>> {
        if !input.starts_with(Self::PREFIX) {
            return Ok(None);
        }

        // Consume the /i segment
        let (i_segment, mut remainder) = input.split_once('/').unwrap_or((input, ""));

        // Check if the next segment is the hydrogen isotope sublayer (/hD2, /hT, etc.)
        // It starts with 'h' followed by D, T, or H (not a digit, which would be the
        // main hydrogen layer pattern like "h1H2").
        let h_isotope_segment = if remainder.starts_with('h')
            && remainder.as_bytes().get(1).is_some_and(|&b| b == b'D' || b == b'T' || b == b'H')
        {
            let (h_seg, rest) = remainder.split_once('/').unwrap_or((remainder, ""));
            remainder = rest;
            Some(h_seg)
        } else {
            None
        };

        *input = remainder;
        let (formula, _) = context;
        let layer = Self::from_str_with_context(i_segment, (formula, h_isotope_segment))?;
        Ok(Some(layer))
    }
}

#[cfg(test)]
mod tests {
    use core::str::FromStr;

    use elements_rs::isotopes::HydrogenIsotope;
    use molecular_formulas::InChIFormula;

    use crate::{
        errors::Error, inchi::isotope_layer::IsotopeLayer, traits::parse::FromStrWithContext,
    };

    fn formula(s: &str) -> InChIFormula {
        InChIFormula::from_str(s).expect("valid formula")
    }

    fn parse(i_layer: &str, h_seg: Option<&str>, ctx: &str) -> Result<IsotopeLayer, Error<u16>> {
        IsotopeLayer::from_str_with_context(i_layer, (&formula(ctx), h_seg))
    }

    #[test]
    fn test_deuterated_water() {
        // i/hD2 → body empty, h_segment = "hD2"
        let result = parse("i", Some("hD2"), "H2O").unwrap();
        assert_eq!(result.components.len(), 1);
        assert!(result.components[0].atoms.is_empty());
        assert_eq!(result.components[0].hydrogens.len(), 1);
        assert_eq!(result.components[0].hydrogens[0].isotope, HydrogenIsotope::D);
        assert_eq!(result.components[0].hydrogens[0].count, 2);
    }

    #[test]
    fn test_tritium() {
        let result = parse("i", Some("hT"), "H2O").unwrap();
        assert_eq!(result.components[0].hydrogens.len(), 1);
        assert_eq!(result.components[0].hydrogens[0].isotope, HydrogenIsotope::T);
        assert_eq!(result.components[0].hydrogens[0].count, 1);
    }

    #[test]
    fn test_protium() {
        let result = parse("i", Some("hH"), "H2O").unwrap();
        assert_eq!(result.components[0].hydrogens.len(), 1);
        assert_eq!(result.components[0].hydrogens[0].isotope, HydrogenIsotope::H1);
        assert_eq!(result.components[0].hydrogens[0].count, 1);
    }

    #[test]
    fn test_atom_mass_shift_positive() {
        // i1+1 → atom 0, mass shift +1
        let result = parse("i1+1", None, "CH4").unwrap();
        assert_eq!(result.components.len(), 1);
        assert_eq!(result.components[0].atoms.len(), 1);
        assert_eq!(result.components[0].atoms[0].atom_index, 0);
        assert_eq!(result.components[0].atoms[0].mass_shift, 1);
    }

    #[test]
    fn test_two_atom_specs() {
        // i1+1,2-1 → atom 0 shift +1, atom 1 shift -1
        let result = parse("i1+1,2-1", None, "C2H6").unwrap();
        assert_eq!(result.components[0].atoms.len(), 2);
        assert_eq!(result.components[0].atoms[0].atom_index, 0);
        assert_eq!(result.components[0].atoms[0].mass_shift, 1);
        assert_eq!(result.components[0].atoms[1].atom_index, 1);
        assert_eq!(result.components[0].atoms[1].mass_shift, -1);
    }

    #[test]
    fn test_multi_component_first_empty() {
        // i;1+2 → component 0 empty, component 1 atom 0 shift +2
        let result = parse("i;1+2", None, "H2O.CH4").unwrap();
        assert_eq!(result.components.len(), 2);
        assert!(result.components[0].atoms.is_empty());
        assert_eq!(result.components[1].atoms.len(), 1);
        assert_eq!(result.components[1].atoms[0].atom_index, 0);
        assert_eq!(result.components[1].atoms[0].mass_shift, 2);
    }

    #[test]
    fn test_empty_body_no_h() {
        // i with nothing → 1 component, no atoms, no H isotopes
        let result = parse("i", None, "H2O").unwrap();
        assert_eq!(result.components.len(), 1);
        assert!(result.components[0].atoms.is_empty());
        assert!(result.components[0].hydrogens.is_empty());
    }

    #[test]
    fn test_wrong_prefix() {
        let err = parse("q+1", None, "H2O").unwrap_err();
        assert!(matches!(err, Error::WrongPrefix));
    }

    #[test]
    fn test_invalid_isotope_char() {
        let err = parse("iX", None, "CH4").unwrap_err();
        assert!(matches!(err, Error::InvalidIsotopeValue('X')));
    }

    #[test]
    fn test_multiple_h_isotopes() {
        // hDT3 → D×1 then T×3
        let result = parse("i", Some("hDT3"), "H2O").unwrap();
        assert_eq!(result.components[0].hydrogens.len(), 2);
        assert_eq!(result.components[0].hydrogens[0].isotope, HydrogenIsotope::D);
        assert_eq!(result.components[0].hydrogens[0].count, 1);
        assert_eq!(result.components[0].hydrogens[1].isotope, HydrogenIsotope::T);
        assert_eq!(result.components[0].hydrogens[1].count, 3);
    }

    #[test]
    fn test_zero_atom_index() {
        // 0 is a leading zero → rejected by try_fold_number
        let err = parse("i0+1", None, "CH4").unwrap_err();
        assert!(matches!(err, Error::InvalidIsotopeValue('0')));
    }

    #[test]
    fn test_missing_sign_in_atom_spec() {
        // "i1" has no +/- sign → error
        let err = parse("i1", None, "CH4").unwrap_err();
        assert!(matches!(err, Error::InvalidIsotopeValue(_)));
    }

    #[test]
    fn test_atom_specs_and_h_isotopes_together() {
        // Atom mass shift + hydrogen isotope sublayer
        let result = parse("i1+1", Some("hD2"), "CH4").unwrap();
        assert_eq!(result.components.len(), 1);
        assert_eq!(result.components[0].atoms.len(), 1);
        assert_eq!(result.components[0].atoms[0].atom_index, 0);
        assert_eq!(result.components[0].atoms[0].mass_shift, 1);
        assert_eq!(result.components[0].hydrogens.len(), 1);
        assert_eq!(result.components[0].hydrogens[0].isotope, HydrogenIsotope::D);
        assert_eq!(result.components[0].hydrogens[0].count, 2);
    }

    #[test]
    fn test_invalid_char_in_h_segment() {
        let err = parse("i", Some("hX"), "H2O").unwrap_err();
        assert!(matches!(err, Error::InvalidIsotopeValue('X')));
    }

    #[test]
    fn test_component_count_mismatch_too_many() {
        // Formula has 1 component but isotope layer has 2
        let err = parse("i1+1;2+1", None, "CH4").unwrap_err();
        assert!(matches!(err, Error::FormulaAndConnectionLayerMixtureMismatch(_)));
    }

    #[test]
    fn test_component_count_mismatch_too_few() {
        // Formula has 2 components but isotope layer has 1
        let err = parse("i1+1", None, "CH4.C2H6").unwrap_err();
        assert!(matches!(err, Error::FormulaAndConnectionLayerMixtureMismatch(_)));
    }

    #[test]
    fn test_repetition_prefix() {
        // 2*1+1 with a 2-component formula → both have atom 0 shift +1
        let result = parse("i2*1+1", None, "2CH4").unwrap();
        assert_eq!(result.components.len(), 2);
        for comp in &result.components {
            assert_eq!(comp.atoms.len(), 1);
            assert_eq!(comp.atoms[0].atom_index, 0);
            assert_eq!(comp.atoms[0].mass_shift, 1);
        }
    }

    #[test]
    fn test_negative_mass_shift() {
        let result = parse("i1-3", None, "CH4").unwrap();
        assert_eq!(result.components[0].atoms[0].mass_shift, -3);
    }

    #[test]
    fn test_large_mass_shift() {
        let result = parse("i1+100", None, "CH4").unwrap();
        assert_eq!(result.components[0].atoms[0].mass_shift, 100);
    }

    #[test]
    fn test_try_build_layer_with_h_sublayer() {
        // Simulates "i1+1/hD2/f..." — should consume both i and h segments
        let f = formula("CH4");
        let mut input = "i1+1/hD2/f1";
        let result = IsotopeLayer::try_build_layer(&mut input, (&f, None)).unwrap().unwrap();
        assert_eq!(input, "f1"); // only /f remains
        assert_eq!(result.components[0].atoms.len(), 1);
        assert_eq!(result.components[0].hydrogens.len(), 1);
        assert_eq!(result.components[0].hydrogens[0].isotope, HydrogenIsotope::D);
    }

    #[test]
    fn test_try_build_layer_does_not_consume_main_h_layer() {
        // "i/h1H2" — the /h segment starts with a digit, so it's the main H layer
        let f = formula("H2O");
        let mut input = "i/h1H2";
        let result = IsotopeLayer::try_build_layer(&mut input, (&f, None)).unwrap().unwrap();
        assert_eq!(input, "h1H2"); // /h1H2 is NOT consumed
        assert!(result.components[0].hydrogens.is_empty());
    }

    #[test]
    fn test_try_build_layer_not_present() {
        let f = formula("CH4");
        let mut input = "f1";
        let result = IsotopeLayer::try_build_layer(&mut input, (&f, None)).unwrap();
        assert!(result.is_none());
        assert_eq!(input, "f1"); // input unchanged
    }

    #[test]
    fn test_try_build_layer_h_segment_at_end() {
        // "i/hD2" with no trailing segments
        let f = formula("H2O");
        let mut input = "i/hD2";
        let result = IsotopeLayer::try_build_layer(&mut input, (&f, None)).unwrap().unwrap();
        assert_eq!(input, ""); // fully consumed
        assert_eq!(result.components[0].hydrogens[0].isotope, HydrogenIsotope::D);
        assert_eq!(result.components[0].hydrogens[0].count, 2);
    }

    // --- h isotope segment edge cases ---

    #[test]
    fn test_h_segment_empty_body() {
        // "h" with no isotope letters → empty hydrogen list
        let result = parse("i", Some("h"), "H2O").unwrap();
        assert!(result.components[0].hydrogens.is_empty());
    }

    #[test]
    fn test_h_segment_leading_zero_count() {
        // "hD02" → leading zero in count is rejected
        let err = parse("i", Some("hD02"), "H2O").unwrap_err();
        assert!(matches!(err, Error::InvalidIsotopeValue('0')));
    }

    #[test]
    fn test_h_segment_trailing_garbage() {
        // "hD2X" → D×2 parsed, then 'X' hits the match
        let err = parse("i", Some("hD2X"), "H2O").unwrap_err();
        assert!(matches!(err, Error::InvalidIsotopeValue('X')));
    }

    // --- atom spec edge cases ---

    #[test]
    fn test_atom_spec_trailing_comma() {
        // "i1+1," → trailing comma produces an empty spec
        let err = parse("i1+1,", None, "CH4").unwrap_err();
        assert!(matches!(err, Error::InvalidIsotopeValue(_)));
    }

    #[test]
    fn test_atom_spec_sign_without_magnitude() {
        // "i1+" → sign but no digits after it → delegated to parse_charge
        let err = parse("i1+", None, "CH4").unwrap_err();
        assert!(matches!(err, Error::InvalidChargeValue(_)));
    }

    #[test]
    fn test_atom_spec_leading_zero_index() {
        // "i01+1" → leading zero before a digit (01 is not valid)
        let err = parse("i01+1", None, "CH4").unwrap_err();
        assert!(matches!(err, Error::InvalidIsotopeValue('0')));
    }

    // --- empty body + multi-component ---

    #[test]
    fn test_empty_body_multi_component_single_result() {
        // Empty atom body with multi-component formula → shortcut returns 1 component.
        // This is the documented behavior: empty body + h sublayer is a
        // single-component shorthand. Multi-component isotope layers use
        // semicolons in the atom body.
        let result = parse("i", Some("hD2"), "H2O.CH4").unwrap();
        assert_eq!(result.components.len(), 1);
    }

    // --- try_build_layer /h disambiguation edge cases ---

    #[test]
    fn test_try_build_layer_bare_h_not_consumed() {
        // "i/h" — 'h' at end with nothing after it; get(1) returns None
        let f = formula("H2O");
        let mut input = "i/h";
        let result = IsotopeLayer::try_build_layer(&mut input, (&f, None)).unwrap().unwrap();
        assert_eq!(input, "h"); // bare /h is NOT consumed as isotope sublayer
        assert!(result.components[0].hydrogens.is_empty());
    }

    #[test]
    fn test_try_build_layer_h_with_slash_not_consumed() {
        // "i/h/f1" — 'h' followed by '/', get(1) is b'/' → not D/T/H
        let f = formula("H2O");
        let mut input = "i/h/f1";
        let result = IsotopeLayer::try_build_layer(&mut input, (&f, None)).unwrap().unwrap();
        assert_eq!(input, "h/f1"); // /h/f1 is NOT consumed
        assert!(result.components[0].hydrogens.is_empty());
    }

    #[test]
    fn test_try_build_layer_i_at_end_of_input() {
        // Just "i" with nothing after — no slash at all
        let f = formula("H2O");
        let mut input = "i";
        let result = IsotopeLayer::try_build_layer(&mut input, (&f, None)).unwrap().unwrap();
        assert_eq!(input, "");
        assert!(result.components[0].atoms.is_empty());
        assert!(result.components[0].hydrogens.is_empty());
    }

    use crate::traits::parse::PrefixFromStrWithContext;
}
