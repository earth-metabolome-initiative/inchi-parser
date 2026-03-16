use molecular_formulas::{BaselineDigit, InChIFormula, MolecularFormula, try_fold_number};

use crate::{
    errors::Error,
    inchi::charge_layer::ChargeSubLayer,
    traits::{
        parse::{FromStrWithContext, PrefixFromStrWithContext},
        prefix::Prefix,
    },
};

/// Parses a charge value from a character iterator.
///
/// Expected formats: empty → 0, `+N` → +N, `-N` → −N.
/// InChI charges always carry an explicit sign for nonzero values.
pub(crate) fn parse_charge(s: &str) -> Result<i16, Error<u16>> {
    if s.is_empty() {
        return Ok(0);
    }

    let mut chars = s.chars().peekable();

    let negative = match chars.peek() {
        Some('+') => {
            chars.next();
            false
        }
        Some('-') => {
            chars.next();
            true
        }
        Some(&c) => return Err(Error::InvalidChargeValue(c)),
        None => unreachable!(),
    };

    let Some(Ok(magnitude)) = try_fold_number::<u16, BaselineDigit, _>(&mut chars) else {
        let c = chars.peek().copied().unwrap_or(if negative { '-' } else { '+' });
        return Err(Error::InvalidChargeValue(c));
    };

    if let Some(&c) = chars.peek() {
        return Err(Error::InvalidChargeValue(c));
    }

    if negative {
        i16::try_from(magnitude).map(|v| -v).or({
            // magnitude == 32768 is valid as -32768 (i16::MIN)
            if magnitude == 32768 { Ok(i16::MIN) } else { Err(Error::InvalidChargeValue('-')) }
        })
    } else {
        i16::try_from(magnitude).map_err(|_| Error::InvalidChargeValue('+'))
    }
}

impl FromStrWithContext for ChargeSubLayer {
    type Context<'a> = &'a InChIFormula;
    type Input<'a> = &'a str;
    type Idx = u16;

    fn from_str_with_context(
        input: Self::Input<'_>,
        context: Self::Context<'_>,
    ) -> Result<Self, Error<Self::Idx>> {
        let s = input.strip_prefix(Self::PREFIX).ok_or(Error::WrongPrefix)?;

        let mut subformulas = context.subformulas();
        let mut charges = alloc::vec::Vec::with_capacity(context.number_of_mixtures());

        for component_str in s.split(';') {
            // Detect optional n* repetition prefix (e.g. "2*+1" → reps=2, rest="+1").
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

            // Consume one subformula for the first occurrence.
            subformulas.next().ok_or(Error::FormulaAndConnectionLayerMixtureMismatch(
                context.number_of_mixtures(),
            ))?;

            let charge = parse_charge(component_str)?;

            charges.push(charge);
            for _ in 1..reps {
                subformulas.next().ok_or(Error::FormulaAndConnectionLayerMixtureMismatch(
                    context.number_of_mixtures(),
                ))?;
                charges.push(charge);
            }
        }

        if subformulas.next().is_some() {
            return Err(Error::FormulaAndConnectionLayerMixtureMismatch(
                context.number_of_mixtures(),
            ));
        }

        Ok(ChargeSubLayer { charges })
    }
}

impl PrefixFromStrWithContext for ChargeSubLayer {}

#[cfg(test)]
mod tests {
    use core::str::FromStr;

    use molecular_formulas::InChIFormula;

    use crate::{
        errors::Error, inchi::charge_layer::ChargeSubLayer, traits::parse::FromStrWithContext,
    };

    fn formula(s: &str) -> InChIFormula {
        InChIFormula::from_str(s).expect("valid formula")
    }

    fn parse(q_layer: &str, ctx: &str) -> Result<ChargeSubLayer, Error<u16>> {
        ChargeSubLayer::from_str_with_context(q_layer, &formula(ctx))
    }

    #[test]
    fn test_single_negative_charge() {
        let result = parse("q-2", "O").unwrap();
        assert_eq!(result.charges, &[-2]);
    }

    #[test]
    fn test_single_positive_charge() {
        let result = parse("q+1", "Na").unwrap();
        assert_eq!(result.charges, &[1]);
    }

    #[test]
    fn test_two_components_neutral_and_positive() {
        // ClH.Na → component 0 neutral, component 1 charge +1
        let result = parse("q;+1", "ClH.Na").unwrap();
        assert_eq!(result.charges, &[0, 1]);
    }

    #[test]
    fn test_two_components_nickel_complex() {
        // Something.Ni → second component charge +2
        let result = parse("q;+2", "O.Ni").unwrap();
        assert_eq!(result.charges, &[0, 2]);
    }

    #[test]
    fn test_repetition_prefix() {
        // 2*+1 with a 2-component formula → [+1, +1]
        let result = parse("q2*+1", "2Na").unwrap();
        assert_eq!(result.charges, &[1, 1]);
    }

    #[test]
    fn test_empty_body_implies_zero() {
        // "q" with nothing after prefix → single empty component → charge 0
        let result = parse("q", "O").unwrap();
        assert_eq!(result.charges, &[0]);
    }

    #[test]
    fn test_mixed_positive_and_negative() {
        let result = parse("q+1;-2", "Na.O").unwrap();
        assert_eq!(result.charges, &[1, -2]);
    }

    #[test]
    fn test_repetition_prefix_negative() {
        let result = parse("q2*-1", "2O").unwrap();
        assert_eq!(result.charges, &[-1, -1]);
    }

    #[test]
    fn test_wrong_prefix() {
        let err = parse("h-2", "O").unwrap_err();
        assert!(matches!(err, Error::WrongPrefix));
    }

    #[test]
    fn test_invalid_charge_value() {
        let err = parse("qabc", "O").unwrap_err();
        assert!(matches!(err, Error::InvalidChargeValue('a')));
    }

    #[test]
    fn test_sign_without_digits() {
        let err = parse("q+", "O").unwrap_err();
        assert!(matches!(err, Error::InvalidChargeValue(_)));
    }

    #[test]
    fn test_component_count_mismatch_too_many() {
        // Formula has 1 component but charge layer has 2
        let err = parse("q+1;+2", "O").unwrap_err();
        assert!(matches!(err, Error::FormulaAndConnectionLayerMixtureMismatch(_)));
    }

    #[test]
    fn test_component_count_mismatch_too_few() {
        // Formula has 2 components but charge layer has 1
        let err = parse("q+1", "O.Na").unwrap_err();
        assert!(matches!(err, Error::FormulaAndConnectionLayerMixtureMismatch(_)));
    }
}
