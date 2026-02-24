use core::str::Chars;

use crate::{
    inchi::main_layer::{AtomConnectionLayer, MolecularGraph},
    traits::{
        parse::{FromStrWithContext, PrefixFromStrWithContext},
        prefix::Prefix,
    },
};
pub mod connection_layer_token_iter;
mod from_connection_layer_token;
use alloc::{format, vec::Vec};

use from_connection_layer_token::FromConnectionLayer;
use geometric_traits::prelude::*;
use molecular_formulas::{BaselineDigit, InChIFormula, MolecularFormula, try_fold_number};

use crate::errors::Error;

impl FromStrWithContext for AtomConnectionLayer {
    type Context<'a> = &'a InChIFormula;
    type Input<'a> = &'a str;
    fn from_str_with_context(
        input: Self::Input<'_>,
        context: Self::Context<'_>,
    ) -> Result<Self, Error<usize>> {
        let Some(s) = input.strip_prefix(Self::PREFIX) else {
            return Err(Error::MissingInchiPrefix);
        };

        // If there are multiple molecules in the molecular formula, we need to split
        // the atom connections layer at the ';' character
        let mut subformulas = context.subformulas();
        let mut molecular_graphs = Vec::with_capacity(context.number_of_mixtures());
        for molecular_graph_input in s.split(';') {
            let mut molecular_graph_chars = molecular_graph_input.chars().peekable();

            let number_of_repetitions = if let Some(Ok(repetitions)) =
                try_fold_number::<u32, BaselineDigit, _>(&mut molecular_graph_chars)
                && molecular_graph_chars.next() == Some('*')
            {
                repetitions
            } else {
                molecular_graph_chars = molecular_graph_input.chars().peekable();
                1
            };

            let Some(subformula) = subformulas.next() else {
                todo!();
                // return Err(Error::FormulaAndConnectionLayerMixtureMismatch(, ));
            };

            let mg = MolecularGraph::from_str_with_context(molecular_graph_chars, &subformula)?;
            for _ in 0..number_of_repetitions {
                let Some(subformula_successor) = subformulas.next() else { todo!() };
                if subformula != subformula_successor {
                    todo!()
                }

                molecular_graphs.push(mg.clone());
            }
        }

        Ok(molecular_graphs) // TODO: if the number of repetition is too big, then this might break ?
    }
}

impl FromStrWithContext for MolecularGraph<usize> {
    type Context<'a> = &'a InChIFormula;
    type Input<'a> = core::iter::Peekable<Chars<'a>>;
    fn from_str_with_context(
        input: Self::Input<'_>,
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
