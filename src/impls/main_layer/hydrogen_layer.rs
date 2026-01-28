use molecular_formulas::InChIFormula;
pub(crate) mod tokens;

use crate::{
    errors::Error,
    inchi::main_layer::HydrogensSubLayer,
    traits::parse::{FromStrWithContext, PrefixFromStrWithContext},
};

impl FromStrWithContext for HydrogensSubLayer {
    type Context<'a> = &'a InChIFormula;
    fn from_str_with_context(
        input: &str,
        context: Self::Context<'_>,
    ) -> Result<Self, Error<usize>> {
        Err(Error::UnimplementedFeature("HydrogenSublayer not implemented yet"))
    }
}

impl PrefixFromStrWithContext for HydrogensSubLayer {}
