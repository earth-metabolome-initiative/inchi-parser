pub mod connection_layer_base_token_iter;
use crate::inchi::main_layer::{AtomConnectionLayer, MolecularGraph};
use crate::traits::parse::ParseLayer;
use crate::traits::prefix::Prefix;
mod from_connection_layer_token;
use from_connection_layer_token::FromConnectionLayer;
use geometric_traits::prelude::*;
use molecular_formulas::MolecularFormula;
use std::str::FromStr;

impl ParseLayer for MolecularFormula {
    type Context<'a> = ();
    type Error = crate::errors::Error;
    fn parse(input: &str, _context: Self::Context<'_>) -> Result<Self, Self::Error> {
        // parse the molecular formula
        let molecular_formula = MolecularFormula::from_str(input)?;
        let is_hill_sorted = molecular_formula.is_hill_sorted()?;
        if !is_hill_sorted {
            return Err(Self::Error::NotHillSorted);
        }
        Ok(molecular_formula)
    }
}

impl ParseLayer for AtomConnectionLayer {
    type Context<'a> = &'a MolecularFormula;
    type Error = crate::errors::Error;
    fn parse(input: &str, context: Self::Context<'_>) -> Result<Self, Self::Error> {
        let Some(s) = input.strip_prefix(Self::PREFIX) else {
            return Ok(None);
        };

        // Check if the connection layer and the molecular formula have the same number of mixtures
        let number_of_mixtures = context.number_of_mixtures();
        let connections_parts = s.matches(";").count() + 1;
        if number_of_mixtures != connections_parts {
            return Err(Self::Error::FormulaAndConnectionLayerMixtureMismatch(
                number_of_mixtures,
                connections_parts,
            ));
        }

        // If there are multiple molecules in the molecular formula, we need to split the atom connections layer
        // at the ';' character

        Ok(Some(
            s.split(';')
                .zip(context.mixtures())
                .map(|(input, mixture)| MolecularGraph::parse(input, mixture))
                .collect::<Result<Vec<_>, _>>()?,
        ))
    }
}

impl ParseLayer for MolecularGraph {
    type Context<'a> = &'a MolecularFormula;
    type Error = crate::errors::Error;
    fn parse(input: &str, context: Self::Context<'_>) -> Result<Self, Self::Error> {
        let vocab_size = context.number_of_elements()?;
        // TODO! REMOVE HYDROGENS!
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
