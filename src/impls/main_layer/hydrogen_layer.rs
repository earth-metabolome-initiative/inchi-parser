use molecular_formulas::{InChIFormula, MolecularFormula};
pub(crate) mod sub_tokens;
mod token_iter;

use crate::{
    errors::Error,
    inchi::main_layer::HydrogensSubLayer,
    traits::{
        parse::{FromStrWithContext, PrefixFromStrWithContext},
        prefix::Prefix,
    },
};

use token_iter::parse_component;

impl FromStrWithContext for HydrogensSubLayer {
    type Context<'a> = &'a InChIFormula;
    type Input<'a> = &'a str;
    type Idx = u16;

    fn from_str_with_context(
        input: Self::Input<'_>,
        context: Self::Context<'_>,
    ) -> Result<Self, Error<Self::Idx>> {
        let s = input
            .strip_prefix(Self::PREFIX)
            .ok_or(Error::WrongPrefix)?;

        let mut subformulas = context.subformulas();
        let mut components = alloc::vec::Vec::with_capacity(context.number_of_mixtures());

        for component_str in s.split(';') {
            // Detect optional n* repetition prefix (e.g. "2*1H" → reps=2, rest="1H").
            // The prefix is purely ASCII digits followed by '*', so we can scan bytes.
            let digit_count = component_str
                .bytes()
                .take_while(u8::is_ascii_digit)
                .count();
            let (reps, component_str) =
                if digit_count > 0 && component_str.as_bytes().get(digit_count) == Some(&b'*') {
                    component_str[..digit_count].parse::<u32>().ok().map_or(
                        (1u32, component_str),
                        |n| (n, &component_str[digit_count + 1..]),
                    )
                } else {
                    (1u32, component_str)
                };

            // Consume one subformula from context for the first occurrence
            let sf = subformulas
                .next()
                .ok_or(Error::FormulaAndConnectionLayerMixtureMismatch(
                    context.number_of_mixtures(),
                ))?;
            let num_atoms = sf.number_of_non_hydrogens();

            // Parse the component
            let component = parse_component::<u16>(component_str, num_atoms)?;

            // Push the first occurrence then n-1 clones for repeated fragments
            components.push(component.clone());
            for _ in 1..reps {
                subformulas
                    .next()
                    .ok_or(Error::FormulaAndConnectionLayerMixtureMismatch(
                        context.number_of_mixtures(),
                    ))?;
                components.push(component.clone());
            }
        }

        if subformulas.next().is_some() {
            return Err(Error::FormulaAndConnectionLayerMixtureMismatch(
                context.number_of_mixtures(),
            ));
        }

        Ok(HydrogensSubLayer { components })
    }
}

impl PrefixFromStrWithContext for HydrogensSubLayer {}

#[cfg(test)]
mod tests {
    use core::str::FromStr;

    use molecular_formulas::InChIFormula;

    use crate::{
        errors::{Error, HydrogenLayerTokenError},
        inchi::main_layer::HydrogensSubLayer,
        traits::parse::FromStrWithContext,
    };

    fn formula(s: &str) -> InChIFormula {
        InChIFormula::from_str(s).expect("valid formula")
    }

    fn parse(h_layer: &str, ctx: &str) -> Result<HydrogensSubLayer, Error<u16>> {
        HydrogensSubLayer::from_str_with_context(h_layer, &formula(ctx))
    }

    #[test]
    fn test_single_fixed_atom() {
        // h1H2 → atom 1 (0-based index 0) has 2 fixed H
        let result = parse("h1H2", "CH4").unwrap();
        assert_eq!(result.components.len(), 1);
        assert_eq!(result.components[0].fixed_h[0], 2);
        assert!(result.components[0].mobile_groups.is_empty());
    }

    #[test]
    fn test_fixed_list_and_range() {
        // h3H,1-2H3 → atom 3 has 1H; atoms 1-2 each have 3H
        let result = parse("h3H,1-2H3", "C3H10").unwrap();
        let fh = &result.components[0].fixed_h;
        assert_eq!(fh[0], 3); // atom 1 → 0-based 0
        assert_eq!(fh[1], 3); // atom 2 → 0-based 1
        assert_eq!(fh[2], 1); // atom 3 → 0-based 2
    }

    #[test]
    fn test_range_all_fixed() {
        // h7-10H → atoms 7-10 each have 1 fixed H
        let result = parse("h7-10H", "C10H22").unwrap();
        let fh = &result.components[0].fixed_h;
        for h in &fh[..6] {
            assert_eq!(*h, 0);
        }
        for h in &fh[6..10] {
            assert_eq!(*h, 1);
        }
    }

    #[test]
    fn test_mobile_then_fixed() {
        // h(H,1,3)2H → mobile {count:1, atoms:[0,2]}; atom 2 (0-based 1) has 1H
        let result = parse("h(H,1,3)2H", "C3H8").unwrap();
        let comp = &result.components[0];
        assert_eq!(comp.fixed_h[1], 1);
        assert_eq!(comp.mobile_groups.len(), 1);
        assert_eq!(comp.mobile_groups[0].count, 1);
        assert!(!comp.mobile_groups[0].charged);
        assert_eq!(&comp.mobile_groups[0].atoms[..], &[0u16, 2u16]);
    }

    #[test]
    fn test_large_mobile_group() {
        // h(H4,10,11,13,14,17) → mobile {count:4, atoms:[9,10,12,13,16]}
        let result = parse("h(H4,10,11,13,14,17)", "C17H36").unwrap();
        let comp = &result.components[0];
        assert!(comp.fixed_h.iter().all(|&h| h == 0));
        assert_eq!(comp.mobile_groups.len(), 1);
        let mg = &comp.mobile_groups[0];
        assert_eq!(mg.count, 4);
        assert!(!mg.charged);
        assert_eq!(&mg.atoms[..], &[9u16, 10u16, 12u16, 13u16, 16u16]);
    }

    #[test]
    fn test_two_adjacent_mobile_groups() {
        // h(H2,1,3)(H2,2,4) → two mobile groups, no comma between them
        let result = parse("h(H2,1,3)(H2,2,4)", "C4H10").unwrap();
        let comp = &result.components[0];
        assert_eq!(comp.mobile_groups.len(), 2);
        assert_eq!(comp.mobile_groups[0].count, 2);
        assert_eq!(&comp.mobile_groups[0].atoms[..], &[0u16, 2u16]);
        assert_eq!(comp.mobile_groups[1].count, 2);
        assert_eq!(&comp.mobile_groups[1].atoms[..], &[1u16, 3u16]);
    }

    #[test]
    fn test_charged_mobile_group() {
        // h(H-,1,2) → mobile {count:1, charged:true, atoms:[0,1]}
        let result = parse("h(H-,1,2)", "C2H6").unwrap();
        let comp = &result.components[0];
        assert_eq!(comp.mobile_groups.len(), 1);
        let mg = &comp.mobile_groups[0];
        assert_eq!(mg.count, 1);
        assert!(mg.charged);
        assert_eq!(&mg.atoms[..], &[0u16, 1u16]);
    }

    #[test]
    fn test_two_components_second_empty() {
        // h1H; → component[0] has fixed_h[0]=1; component[1] is empty
        let result = parse("h1H;", "CH4.C2H6").unwrap();
        assert_eq!(result.components.len(), 2);
        assert_eq!(result.components[0].fixed_h[0], 1);
        assert!(result.components[1].fixed_h.iter().all(|&h| h == 0));
        assert!(result.components[1].mobile_groups.is_empty());
    }

    #[test]
    fn test_repetition_prefix() {
        // h2*1H → 2 components each with fixed_h[0] = 1
        let result = parse("h2*1H", "2CH4").unwrap();
        assert_eq!(result.components.len(), 2);
        for comp in &result.components {
            assert_eq!(comp.fixed_h[0], 1);
            assert!(comp.mobile_groups.is_empty());
        }
    }

    #[test]
    fn test_wrong_prefix() {
        let err = parse("c1-2-3", "C3H8").unwrap_err();
        assert!(matches!(err, Error::WrongPrefix));
    }

    #[test]
    fn test_invalid_character() {
        let err = parse("h1X", "CH4").unwrap_err();
        assert!(matches!(
            err,
            Error::HydrogenLayerTokenError(HydrogenLayerTokenError::InvalidCharacter('X'))
        ));
    }

}
