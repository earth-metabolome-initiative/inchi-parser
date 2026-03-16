use alloc::vec::Vec;

use molecular_formulas::{BaselineDigit, InChIFormula, MolecularFormula, try_fold_number};

use crate::{
    errors::Error,
    inchi::stereochemistry_layer::{
        AlleneSublayer, DoubleBondStereo, DoubleBondSublayer, StereoChemistryInformationSublayer,
        StereoParity, TetrahedralStereo, TetrahedralSublayer,
    },
    traits::{
        parse::{FromStrWithContext, PrefixFromStrWithContext},
        prefix::Prefix,
    },
};

/// Parses a parity character (`+`, `-`, `?`) into a `StereoParity`.
fn parse_parity(c: char) -> Result<StereoParity, Error<u16>> {
    match c {
        '+' => Ok(StereoParity::Plus),
        '-' => Ok(StereoParity::Minus),
        '?' => Ok(StereoParity::Unknown),
        c => Err(Error::InvalidStereoValue(c)),
    }
}

/// Parses a single double bond spec like `12-11+`.
fn parse_double_bond_spec(spec: &str) -> Result<DoubleBondStereo, Error<u16>> {
    let mut chars = spec.chars().peekable();

    // Parse first atom index (1-based)
    let atom1_1based = match try_fold_number::<u16, BaselineDigit, _>(&mut chars) {
        Some(Ok(n)) if n > 0 => n,
        _ => {
            return Err(Error::InvalidStereoValue(spec.chars().next().unwrap_or('?')));
        }
    };

    // Consume '-' separator
    match chars.next() {
        Some('-') => {}
        Some(c) => return Err(Error::InvalidStereoValue(c)),
        None => return Err(Error::InvalidStereoValue('-')),
    }

    // Parse second atom index (1-based)
    let atom2_1based = match try_fold_number::<u16, BaselineDigit, _>(&mut chars) {
        Some(Ok(n)) if n > 0 => n,
        _ => {
            return Err(Error::InvalidStereoValue(chars.next().unwrap_or('?')));
        }
    };

    // Parse parity
    let parity = match chars.next() {
        Some(c) => parse_parity(c)?,
        None => return Err(Error::InvalidStereoValue('?')),
    };

    if chars.next().is_some() {
        return Err(Error::InvalidStereoValue('?'));
    }

    Ok(DoubleBondStereo { atom1: atom1_1based - 1, atom2: atom2_1based - 1, parity })
}

/// Parses a single tetrahedral spec like `13-`.
fn parse_tetrahedral_spec(spec: &str) -> Result<TetrahedralStereo, Error<u16>> {
    let mut chars = spec.chars().peekable();

    // Parse atom index (1-based)
    let atom_1based = match try_fold_number::<u16, BaselineDigit, _>(&mut chars) {
        Some(Ok(n)) if n > 0 => n,
        _ => {
            return Err(Error::InvalidStereoValue(spec.chars().next().unwrap_or('?')));
        }
    };

    // Parse parity
    let parity = match chars.next() {
        Some(c) => parse_parity(c)?,
        None => return Err(Error::InvalidStereoValue('?')),
    };

    if chars.next().is_some() {
        return Err(Error::InvalidStereoValue('?'));
    }

    Ok(TetrahedralStereo { atom: atom_1based - 1, parity })
}

// --- DoubleBondSublayer ---

impl FromStrWithContext for DoubleBondSublayer {
    type Context<'a> = &'a InChIFormula;
    type Input<'a> = &'a str;
    type Idx = u16;

    fn from_str_with_context(
        input: Self::Input<'_>,
        context: Self::Context<'_>,
    ) -> Result<Self, Error<Self::Idx>> {
        let s = input.strip_prefix(Self::PREFIX).ok_or(Error::WrongPrefix)?;

        let mut subformulas = context.subformulas();
        let mut components = Vec::with_capacity(context.number_of_mixtures());

        for component_str in s.split(';') {
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
                context.number_of_mixtures(),
            ))?;

            let bonds = if component_str.is_empty() {
                Vec::new()
            } else {
                component_str
                    .split(',')
                    .map(parse_double_bond_spec)
                    .collect::<Result<Vec<_>, _>>()?
            };

            components.push(bonds.clone());
            for _ in 1..reps {
                subformulas.next().ok_or(Error::FormulaAndConnectionLayerMixtureMismatch(
                    context.number_of_mixtures(),
                ))?;
                components.push(bonds.clone());
            }
        }

        if subformulas.next().is_some() {
            return Err(Error::FormulaAndConnectionLayerMixtureMismatch(
                context.number_of_mixtures(),
            ));
        }

        Ok(DoubleBondSublayer { components })
    }
}

impl PrefixFromStrWithContext for DoubleBondSublayer {}

// --- TetrahedralSublayer ---

impl FromStrWithContext for TetrahedralSublayer {
    type Context<'a> = &'a InChIFormula;
    type Input<'a> = &'a str;
    type Idx = u16;

    fn from_str_with_context(
        input: Self::Input<'_>,
        context: Self::Context<'_>,
    ) -> Result<Self, Error<Self::Idx>> {
        let s = input.strip_prefix(Self::PREFIX).ok_or(Error::WrongPrefix)?;

        let mut subformulas = context.subformulas();
        let mut components = Vec::with_capacity(context.number_of_mixtures());

        for component_str in s.split(';') {
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
                context.number_of_mixtures(),
            ))?;

            let centers = if component_str.is_empty() {
                Vec::new()
            } else {
                component_str
                    .split(',')
                    .map(parse_tetrahedral_spec)
                    .collect::<Result<Vec<_>, _>>()?
            };

            components.push(centers.clone());
            for _ in 1..reps {
                subformulas.next().ok_or(Error::FormulaAndConnectionLayerMixtureMismatch(
                    context.number_of_mixtures(),
                ))?;
                components.push(centers.clone());
            }
        }

        if subformulas.next().is_some() {
            return Err(Error::FormulaAndConnectionLayerMixtureMismatch(
                context.number_of_mixtures(),
            ));
        }

        Ok(TetrahedralSublayer { components })
    }
}

impl PrefixFromStrWithContext for TetrahedralSublayer {}

// --- AlleneSublayer ---

impl FromStrWithContext for AlleneSublayer {
    type Context<'a> = ();
    type Input<'a> = &'a str;
    type Idx = u16;

    fn from_str_with_context(
        input: Self::Input<'_>,
        _context: Self::Context<'_>,
    ) -> Result<Self, Error<Self::Idx>> {
        let s = input.strip_prefix(Self::PREFIX).ok_or(Error::WrongPrefix)?;

        let mut values = Vec::new();
        for segment in s.split('.') {
            if segment.is_empty() {
                values.push(None);
            } else {
                for c in segment.chars() {
                    match c {
                        '0' => values.push(Some(0)),
                        '1' => values.push(Some(1)),
                        _ => return Err(Error::InvalidStereoValue(c)),
                    }
                }
            }
        }

        Ok(AlleneSublayer { values })
    }
}

impl PrefixFromStrWithContext for AlleneSublayer {}

// --- StereoChemistryInformationSublayer ---

impl FromStrWithContext for StereoChemistryInformationSublayer {
    type Context<'a> = ();
    type Input<'a> = &'a str;
    type Idx = u16;

    fn from_str_with_context(
        input: Self::Input<'_>,
        _context: Self::Context<'_>,
    ) -> Result<Self, Error<Self::Idx>> {
        let s = input.strip_prefix(Self::PREFIX).ok_or(Error::WrongPrefix)?;

        let mut chars = s.chars();
        let c = chars.next().ok_or(Error::InvalidStereoValue('s'))?;
        if chars.next().is_some() {
            return Err(Error::InvalidStereoValue(c));
        }
        match c {
            '1'..='3' => Ok(StereoChemistryInformationSublayer { value: c as u8 - b'0' }),
            _ => Err(Error::InvalidStereoValue(c)),
        }
    }
}

impl PrefixFromStrWithContext for StereoChemistryInformationSublayer {}

#[cfg(test)]
mod tests {
    use core::str::FromStr;

    use molecular_formulas::InChIFormula;

    use crate::{
        errors::Error,
        inchi::stereochemistry_layer::{
            AlleneSublayer, DoubleBondSublayer, StereoChemistryInformationSublayer, StereoParity,
            TetrahedralSublayer,
        },
        traits::parse::FromStrWithContext,
    };

    fn formula(s: &str) -> InChIFormula {
        InChIFormula::from_str(s).expect("valid formula")
    }

    // --- /b tests ---

    #[test]
    fn test_b_single_bond_trans() {
        let f = formula("C12H20");
        let result = DoubleBondSublayer::from_str_with_context("b12-11+", &f).unwrap();
        assert_eq!(result.components.len(), 1);
        assert_eq!(result.components[0].len(), 1);
        assert_eq!(result.components[0][0].atom1, 11);
        assert_eq!(result.components[0][0].atom2, 10);
        assert_eq!(result.components[0][0].parity, StereoParity::Plus);
    }

    #[test]
    fn test_b_two_bonds() {
        let f = formula("C17H30");
        let result = DoubleBondSublayer::from_str_with_context("b12-11+,17-13+", &f).unwrap();
        assert_eq!(result.components[0].len(), 2);
        assert_eq!(result.components[0][0].atom1, 11);
        assert_eq!(result.components[0][0].atom2, 10);
        assert_eq!(result.components[0][0].parity, StereoParity::Plus);
        assert_eq!(result.components[0][1].atom1, 16);
        assert_eq!(result.components[0][1].atom2, 12);
        assert_eq!(result.components[0][1].parity, StereoParity::Plus);
    }

    #[test]
    fn test_b_plus_and_unknown() {
        let f = formula("C10H16");
        let result = DoubleBondSublayer::from_str_with_context("b6-5+,10-7?", &f).unwrap();
        assert_eq!(result.components[0].len(), 2);
        assert_eq!(result.components[0][0].parity, StereoParity::Plus);
        assert_eq!(result.components[0][1].parity, StereoParity::Unknown);
    }

    #[test]
    fn test_b_minus_parity() {
        let f = formula("C11H20");
        let result = DoubleBondSublayer::from_str_with_context("b11-6-", &f).unwrap();
        assert_eq!(result.components[0][0].parity, StereoParity::Minus);
    }

    #[test]
    fn test_b_multi_component() {
        let f = formula("C3H6.C2H4");
        let result = DoubleBondSublayer::from_str_with_context("b3-2+;1-2-", &f).unwrap();
        assert_eq!(result.components.len(), 2);
        assert_eq!(result.components[0].len(), 1);
        assert_eq!(result.components[0][0].atom1, 2);
        assert_eq!(result.components[0][0].atom2, 1);
        assert_eq!(result.components[0][0].parity, StereoParity::Plus);
        assert_eq!(result.components[1].len(), 1);
        assert_eq!(result.components[1][0].atom1, 0);
        assert_eq!(result.components[1][0].atom2, 1);
        assert_eq!(result.components[1][0].parity, StereoParity::Minus);
    }

    #[test]
    fn test_b_empty_body() {
        let f = formula("C3H6");
        let result = DoubleBondSublayer::from_str_with_context("b", &f).unwrap();
        assert_eq!(result.components.len(), 1);
        assert!(result.components[0].is_empty());
    }

    #[test]
    fn test_b_wrong_prefix() {
        let f = formula("C3H6");
        let err = DoubleBondSublayer::from_str_with_context("t3-2+", &f).unwrap_err();
        assert!(matches!(err, Error::WrongPrefix));
    }

    // --- /t tests ---

    #[test]
    fn test_t_four_centers() {
        let f = formula("C16H30");
        let result = TetrahedralSublayer::from_str_with_context("t13-,14-,15+,16-", &f).unwrap();
        assert_eq!(result.components[0].len(), 4);
        assert_eq!(result.components[0][0].atom, 12);
        assert_eq!(result.components[0][0].parity, StereoParity::Minus);
        assert_eq!(result.components[0][1].atom, 13);
        assert_eq!(result.components[0][1].parity, StereoParity::Minus);
        assert_eq!(result.components[0][2].atom, 14);
        assert_eq!(result.components[0][2].parity, StereoParity::Plus);
        assert_eq!(result.components[0][3].atom, 15);
        assert_eq!(result.components[0][3].parity, StereoParity::Minus);
    }

    #[test]
    fn test_t_two_centers() {
        let f = formula("C5H10");
        let result = TetrahedralSublayer::from_str_with_context("t2-,5+", &f).unwrap();
        assert_eq!(result.components[0].len(), 2);
        assert_eq!(result.components[0][0].atom, 1);
        assert_eq!(result.components[0][0].parity, StereoParity::Minus);
        assert_eq!(result.components[0][1].atom, 4);
        assert_eq!(result.components[0][1].parity, StereoParity::Plus);
    }

    #[test]
    fn test_t_unknown_parity() {
        let f = formula("C6H12");
        let result = TetrahedralSublayer::from_str_with_context("t6?", &f).unwrap();
        assert_eq!(result.components[0].len(), 1);
        assert_eq!(result.components[0][0].atom, 5);
        assert_eq!(result.components[0][0].parity, StereoParity::Unknown);
    }

    #[test]
    fn test_t_multi_component() {
        let f = formula("C32H34N4O4.Ni");
        let result = TetrahedralSublayer::from_str_with_context("t15-,19-;", &f).unwrap();
        assert_eq!(result.components.len(), 2);
        assert_eq!(result.components[0].len(), 2);
        assert_eq!(result.components[0][0].atom, 14);
        assert_eq!(result.components[0][0].parity, StereoParity::Minus);
        assert_eq!(result.components[0][1].atom, 18);
        assert_eq!(result.components[0][1].parity, StereoParity::Minus);
        assert!(result.components[1].is_empty());
    }

    #[test]
    fn test_t_repetition_prefix() {
        let f = formula("2CH4");
        let result = TetrahedralSublayer::from_str_with_context("t2*1+", &f).unwrap();
        assert_eq!(result.components.len(), 2);
        for comp in &result.components {
            assert_eq!(comp.len(), 1);
            assert_eq!(comp[0].atom, 0);
            assert_eq!(comp[0].parity, StereoParity::Plus);
        }
    }

    #[test]
    fn test_t_empty_body() {
        let f = formula("CH4");
        let result = TetrahedralSublayer::from_str_with_context("t", &f).unwrap();
        assert_eq!(result.components.len(), 1);
        assert!(result.components[0].is_empty());
    }

    #[test]
    fn test_t_wrong_prefix() {
        let f = formula("CH4");
        let err = TetrahedralSublayer::from_str_with_context("b3+", &f).unwrap_err();
        assert!(matches!(err, Error::WrongPrefix));
    }

    // --- /m tests ---

    #[test]
    fn test_m_zero() {
        let result = AlleneSublayer::from_str_with_context("m0", ()).unwrap();
        assert_eq!(result.values, &[Some(0)]);
    }

    #[test]
    fn test_m_one() {
        let result = AlleneSublayer::from_str_with_context("m1", ()).unwrap();
        assert_eq!(result.values, &[Some(1)]);
    }

    #[test]
    fn test_m_zero_dot() {
        let result = AlleneSublayer::from_str_with_context("m0.", ()).unwrap();
        assert_eq!(result.values, &[Some(0), None]);
    }

    #[test]
    fn test_m_one_dot() {
        let result = AlleneSublayer::from_str_with_context("m1.", ()).unwrap();
        assert_eq!(result.values, &[Some(1), None]);
    }

    #[test]
    fn test_m_two_digits_same_group() {
        // m00. → two components value 0 in first group, empty second group
        let result = AlleneSublayer::from_str_with_context("m00.", ()).unwrap();
        assert_eq!(result.values, &[Some(0), Some(0), None]);
    }

    #[test]
    fn test_m_two_digits_no_dot() {
        // m01 → two components: first 0, second 1
        let result = AlleneSublayer::from_str_with_context("m01", ()).unwrap();
        assert_eq!(result.values, &[Some(0), Some(1)]);
    }

    #[test]
    fn test_m_empty_first_group() {
        // m.11 → first group empty, then two components value 1
        let result = AlleneSublayer::from_str_with_context("m.11", ()).unwrap();
        assert_eq!(result.values, &[None, Some(1), Some(1)]);
    }

    #[test]
    fn test_m_two_digits_trailing_dot() {
        // m11. → two components value 1, empty third group
        let result = AlleneSublayer::from_str_with_context("m11.", ()).unwrap();
        assert_eq!(result.values, &[Some(1), Some(1), None]);
    }

    #[test]
    fn test_m_wrong_prefix() {
        let err = AlleneSublayer::from_str_with_context("s1", ()).unwrap_err();
        assert!(matches!(err, Error::WrongPrefix));
    }

    #[test]
    fn test_m_invalid_value() {
        let err = AlleneSublayer::from_str_with_context("m2", ()).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue('2')));
    }

    // --- /s tests ---

    #[test]
    fn test_s_absolute() {
        let result = StereoChemistryInformationSublayer::from_str_with_context("s1", ()).unwrap();
        assert_eq!(result.value, 1);
    }

    #[test]
    fn test_s_relative() {
        let result = StereoChemistryInformationSublayer::from_str_with_context("s2", ()).unwrap();
        assert_eq!(result.value, 2);
    }

    #[test]
    fn test_s_racemic() {
        let result = StereoChemistryInformationSublayer::from_str_with_context("s3", ()).unwrap();
        assert_eq!(result.value, 3);
    }

    #[test]
    fn test_s_wrong_prefix() {
        let err = StereoChemistryInformationSublayer::from_str_with_context("m1", ()).unwrap_err();
        assert!(matches!(err, Error::WrongPrefix));
    }

    #[test]
    fn test_s_invalid_value() {
        let err = StereoChemistryInformationSublayer::from_str_with_context("s0", ()).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue('0')));
    }

    #[test]
    fn test_s_empty_body() {
        let err = StereoChemistryInformationSublayer::from_str_with_context("s", ()).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue('s')));
    }

    #[test]
    fn test_s_value_4_rejected() {
        let err = StereoChemistryInformationSublayer::from_str_with_context("s4", ()).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue('4')));
    }

    #[test]
    fn test_s_multi_digit_rejected() {
        // "s12" → first char '1' consumed, then '2' remains → error
        let err = StereoChemistryInformationSublayer::from_str_with_context("s12", ()).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue('1')));
    }

    // --- /b error cases ---

    #[test]
    fn test_b_missing_parity() {
        // "b12-11" → no parity char after second atom
        let f = formula("C12H20");
        let err = DoubleBondSublayer::from_str_with_context("b12-11", &f).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue(_)));
    }

    #[test]
    fn test_b_missing_separator() {
        // "b1211+" → no '-' between atoms
        let f = formula("C12H20");
        let err = DoubleBondSublayer::from_str_with_context("b1211+", &f).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue(_)));
    }

    #[test]
    fn test_b_zero_atom_index() {
        // "b0-1+" → atom1 is 0, invalid
        let f = formula("C3H6");
        let err = DoubleBondSublayer::from_str_with_context("b0-1+", &f).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue(_)));
    }

    #[test]
    fn test_b_invalid_parity_char() {
        // "b1-2x" → 'x' is not a valid parity
        let f = formula("C3H6");
        let err = DoubleBondSublayer::from_str_with_context("b1-2x", &f).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue('x')));
    }

    #[test]
    fn test_b_component_mismatch() {
        // Two components in /b but single component formula
        let f = formula("C3H6");
        let err = DoubleBondSublayer::from_str_with_context("b1-2+;1-2-", &f).unwrap_err();
        assert!(matches!(err, Error::FormulaAndConnectionLayerMixtureMismatch(_)));
    }

    #[test]
    fn test_b_trailing_comma() {
        // "b1-2+," → trailing comma produces empty spec
        let f = formula("C3H6");
        let err = DoubleBondSublayer::from_str_with_context("b1-2+,", &f).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue(_)));
    }

    // --- /t error cases ---

    #[test]
    fn test_t_missing_parity() {
        // "t13" → no parity char after atom
        let f = formula("C16H30");
        let err = TetrahedralSublayer::from_str_with_context("t13", &f).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue(_)));
    }

    #[test]
    fn test_t_zero_atom_index() {
        // "t0+" → atom index is 0, invalid
        let f = formula("CH4");
        let err = TetrahedralSublayer::from_str_with_context("t0+", &f).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue(_)));
    }

    #[test]
    fn test_t_invalid_parity_char() {
        // "t1x" → 'x' is not a valid parity
        let f = formula("CH4");
        let err = TetrahedralSublayer::from_str_with_context("t1x", &f).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue('x')));
    }

    #[test]
    fn test_t_trailing_comma() {
        // "t1+," → trailing comma produces empty spec
        let f = formula("CH4");
        let err = TetrahedralSublayer::from_str_with_context("t1+,", &f).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue(_)));
    }

    #[test]
    fn test_t_component_mismatch() {
        // Two components in /t but single component formula
        let f = formula("CH4");
        let err = TetrahedralSublayer::from_str_with_context("t1+;2-", &f).unwrap_err();
        assert!(matches!(err, Error::FormulaAndConnectionLayerMixtureMismatch(_)));
    }

    // --- /m edge cases ---

    #[test]
    fn test_m_multiple_dots() {
        // m0..1 → three groups: first has '0', second empty, third has '1'
        let result = AlleneSublayer::from_str_with_context("m0..1", ()).unwrap();
        assert_eq!(result.values, &[Some(0), None, Some(1)]);
    }

    #[test]
    fn test_m_all_empty() {
        // m.. → two dots = three groups all empty
        let result = AlleneSublayer::from_str_with_context("m..", ()).unwrap();
        assert_eq!(result.values, &[None, None, None]);
    }

    #[test]
    fn test_m_letter_rejected() {
        let err = AlleneSublayer::from_str_with_context("ma", ()).unwrap_err();
        assert!(matches!(err, Error::InvalidStereoValue('a')));
    }

    // --- try_build_layer tests ---

    use crate::traits::parse::PrefixFromStrWithContext;

    #[test]
    fn test_b_try_build_layer_not_present() {
        let f = formula("CH4");
        let mut input = "t1+";
        let result = DoubleBondSublayer::try_build_layer(&mut input, &f).unwrap();
        assert!(result.is_none());
        assert_eq!(input, "t1+"); // unchanged
    }

    #[test]
    fn test_t_try_build_layer_consumes_segment() {
        let f = formula("CH4");
        let mut input = "t1+/m0/s1";
        let result = TetrahedralSublayer::try_build_layer(&mut input, &f).unwrap().unwrap();
        assert_eq!(result.components[0].len(), 1);
        assert_eq!(input, "m0/s1"); // /t consumed
    }

    #[test]
    fn test_m_try_build_layer_consumes_segment() {
        let mut input = "m0/s1";
        let result = AlleneSublayer::try_build_layer(&mut input, ()).unwrap().unwrap();
        assert_eq!(result.values, &[Some(0)]);
        assert_eq!(input, "s1"); // /m consumed
    }

    #[test]
    fn test_s_try_build_layer_at_end() {
        let mut input = "s1";
        let result =
            StereoChemistryInformationSublayer::try_build_layer(&mut input, ()).unwrap().unwrap();
        assert_eq!(result.value, 1);
        assert_eq!(input, ""); // fully consumed
    }

    #[test]
    fn test_b_try_build_layer_at_end() {
        let f = formula("C3H6");
        let mut input = "b1-2+";
        let result = DoubleBondSublayer::try_build_layer(&mut input, &f).unwrap().unwrap();
        assert_eq!(result.components[0].len(), 1);
        assert_eq!(input, ""); // fully consumed
    }
}
