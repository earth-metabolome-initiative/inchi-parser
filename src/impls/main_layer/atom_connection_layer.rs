use crate::{
    inchi::main_layer::{AtomConnectionLayer, MolecularGraph},
    traits::{
        parse::{FromStrWithContext, PrefixFromStrWithContext},
        prefix::Prefix,
    },
};
pub mod connection_layer_base_token_iter;
mod from_connection_layer_token;
use alloc::{format, vec::Vec};

use from_connection_layer_token::FromConnectionLayer;
use geometric_traits::prelude::*;
use molecular_formulas::{InChIFormula, MolecularFormula};

use crate::errors::Error;

impl FromStrWithContext for AtomConnectionLayer {
    type Context<'a> = &'a InChIFormula;
    fn from_str_with_context(
        input: &str,
        context: Self::Context<'_>,
    ) -> Result<Self, Error<usize>> {
        let Some(s) = input.strip_prefix(Self::PREFIX) else {
            return Err(Error::MissingInchiPrefix);
        };

        // Check if the connection layer and the molecular formula have the same number
        // of mixtures
        let number_of_mixtures = context.number_of_mixtures();
        let connections_parts = s.matches(";").count() + 1;
        if number_of_mixtures != connections_parts {
            return Err(Error::FormulaAndConnectionLayerMixtureMismatch(
                number_of_mixtures,
                connections_parts,
            ));
        }

        // If there are multiple molecules in the molecular formula, we need to split
        // the atom connections layer at the ';' character

        Ok(s.split(';')
            .zip(context.subformulas())
            .map(|(input, mixture)| MolecularGraph::from_str_with_context(input, &mixture))
            .collect::<Result<Vec<_>, _>>()?)
    }
}

impl FromStrWithContext for MolecularGraph<usize> {
    type Context<'a> = &'a InChIFormula;
    fn from_str_with_context(
        input: &str,
        context: Self::Context<'_>,
    ) -> Result<Self, Error<usize>> {
        let vocab_size = context.number_of_non_hydrogens();
        let mut edges: Vec<(usize, usize)> = Vec::from_connection_layer_token(input)?;
        edges.sort_unstable();

        let edges: SymmetricCSR2D<CSR2D<usize, usize, usize>> = UndiEdgesBuilder::default()
            .expected_number_of_edges(edges.len())
            .expected_shape(vocab_size)
            .edges(edges.clone().into_iter())
            .build()
            .expect(&format!(
                "Failed to build UndiEdgesBuilder: vocabulary size: {}, edges: {:?}, number of edges: {}",
                vocab_size, edges, edges.len()
            ));

        Ok(MolecularGraph::from((vocab_size, edges)))
    }
}

impl PrefixFromStrWithContext for AtomConnectionLayer {}
