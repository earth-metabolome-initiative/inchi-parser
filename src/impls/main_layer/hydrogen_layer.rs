use molecular_formulas::InChIFormula;
pub(crate) mod sub_tokens;
mod token_iter;

use crate::{
    errors::Error,
    inchi::main_layer::HydrogensSubLayer,
    traits::parse::{FromStrWithContext, PrefixFromStrWithContext},
};

impl FromStrWithContext for HydrogensSubLayer {
    type Context<'a> = &'a InChIFormula;
    type Input<'a> = &'a str;
    fn from_str_with_context(
        input: Self::Input<'_>,
        context: Self::Context<'_>,
    ) -> Result<Self, Error<usize>> {
        Err(Error::UnimplementedFeature("HydrogenSublayer not implemented yet"))
    }
}

impl PrefixFromStrWithContext for HydrogensSubLayer {}
