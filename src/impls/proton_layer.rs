use crate::{
    errors::Error,
    impls::charge_layer::parse_charge,
    inchi::proton_layer::ProtonSublayer,
    traits::{
        parse::{FromStrWithContext, PrefixFromStrWithContext},
        prefix::Prefix,
    },
};

impl FromStrWithContext for ProtonSublayer {
    type Context<'a> = ();
    type Input<'a> = &'a str;
    type Idx = u16;

    fn from_str_with_context(
        input: Self::Input<'_>,
        _context: Self::Context<'_>,
    ) -> Result<Self, Error<Self::Idx>> {
        let s = input.strip_prefix(Self::PREFIX).ok_or(Error::WrongPrefix)?;
        let proton_count = parse_charge(s)?;
        Ok(ProtonSublayer { proton_count })
    }
}

impl PrefixFromStrWithContext for ProtonSublayer {}

#[cfg(test)]
mod tests {
    use crate::{
        errors::Error, inchi::proton_layer::ProtonSublayer, traits::parse::FromStrWithContext,
    };

    fn parse(p_layer: &str) -> Result<ProtonSublayer, Error<u16>> {
        ProtonSublayer::from_str_with_context(p_layer, ())
    }

    #[test]
    fn test_negative_one() {
        let result = parse("p-1").unwrap();
        assert_eq!(result.proton_count, -1);
    }

    #[test]
    fn test_negative_two() {
        let result = parse("p-2").unwrap();
        assert_eq!(result.proton_count, -2);
    }

    #[test]
    fn test_positive_one() {
        let result = parse("p+1").unwrap();
        assert_eq!(result.proton_count, 1);
    }

    #[test]
    fn test_positive_large() {
        let result = parse("p+5").unwrap();
        assert_eq!(result.proton_count, 5);
    }

    #[test]
    fn test_empty_body_implies_zero() {
        let result = parse("p").unwrap();
        assert_eq!(result.proton_count, 0);
    }

    #[test]
    fn test_wrong_prefix() {
        let err = parse("q-1").unwrap_err();
        assert!(matches!(err, Error::WrongPrefix));
    }

    #[test]
    fn test_invalid_value() {
        let err = parse("pabc").unwrap_err();
        assert!(matches!(err, Error::InvalidChargeValue('a')));
    }

    #[test]
    fn test_sign_without_digits() {
        let err = parse("p+").unwrap_err();
        assert!(matches!(err, Error::InvalidChargeValue(_)));
    }
}
