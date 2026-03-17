use core::str::Chars;

use crate::{
    inchi::main_layer::{AtomConnectionComponent, AtomConnectionLayer, MolecularGraph},
    traits::{
        parse::{FromStrWithContext, PrefixFromStrWithContext},
        prefix::Prefix,
    },
};
pub mod connection_layer_token_iter;
mod from_connection_layer_token;
use alloc::{string::String, vec::Vec};

use from_connection_layer_token::FromConnectionLayer;
use geometric_traits::prelude::*;
use molecular_formulas::{BaselineDigit, InChIFormula, MolecularFormula, try_fold_number};

use crate::errors::Error;

impl FromStrWithContext for AtomConnectionLayer {
    type Context<'a> = &'a InChIFormula;
    type Input<'a> = &'a str;
    type Idx = u16;
    fn from_str_with_context(
        input: Self::Input<'_>,
        context: Self::Context<'_>,
    ) -> Result<Self, Error<Self::Idx>> {
        let Some(s) = input.strip_prefix(Self::PREFIX) else {
            return Err(Error::MissingInchiPrefix);
        };

        // If there are multiple molecules in the molecular formula, we need to split
        // the atom connections layer at the ';' character
        let mut subformulas = context.subformulas();
        let mut molecular_graphs = Vec::with_capacity(context.number_of_mixtures());
        for molecular_graph_input in s.split(';') {
            let mut molecular_graph_chars = molecular_graph_input.chars().peekable();

            let (number_of_repetitions, connection_body) = if let Some(Ok(repetitions)) =
                try_fold_number::<u32, BaselineDigit, _>(&mut molecular_graph_chars)
                && molecular_graph_chars.next() == Some('*')
            {
                // Body starts after "N*"
                let star_pos = molecular_graph_input.find('*').unwrap();
                (repetitions, &molecular_graph_input[star_pos + 1..])
            } else {
                molecular_graph_chars = molecular_graph_input.chars().peekable();
                (1, molecular_graph_input)
            };

            let Some(subformula) = subformulas.next() else {
                return Err(Error::FormulaAndConnectionLayerMixtureMismatch(
                    context.number_of_mixtures(),
                ));
            };

            let mg: GenericGraph<u16, SymmetricCSR2D<CSR2D<u16, u16, u16>>> =
                MolecularGraph::from_str_with_context(molecular_graph_chars, &subformula)?;
            molecular_graphs.push(AtomConnectionComponent {
                graph: mg,
                connection_string: String::from(connection_body),
            });

            for _ in 1..number_of_repetitions {
                let Some(subformula) = subformulas.next() else {
                    return Err(Error::FormulaAndConnectionLayerMixtureMismatch(
                        context.number_of_mixtures(),
                    ));
                };

                let chars = connection_body.chars().peekable();
                let mg: GenericGraph<u16, SymmetricCSR2D<CSR2D<u16, u16, u16>>> =
                    MolecularGraph::from_str_with_context(chars, &subformula)?;
                molecular_graphs.push(AtomConnectionComponent {
                    graph: mg,
                    connection_string: String::from(connection_body),
                });
            }
        }

        if subformulas.next().is_some() {
            return Err(Error::FormulaAndConnectionLayerMixtureMismatch(
                context.number_of_mixtures(),
            ));
        }

        Ok(molecular_graphs)
    }
}

impl FromStrWithContext for MolecularGraph<u16> {
    type Context<'a> = &'a InChIFormula;
    type Input<'a> = core::iter::Peekable<Chars<'a>>;
    type Idx = u16;
    fn from_str_with_context(
        input: Self::Input<'_>,
        context: Self::Context<'_>,
    ) -> Result<Self, Error<Self::Idx>> {
        let vocab_size = Self::Idx::try_from(context.number_of_elements())?;
        let mut edges: Vec<(Self::Idx, Self::Idx)> = Vec::from_connection_layer_token(input)?;

        // Validate all atom indices are within the formula's atom count.
        for &(a, b) in &edges {
            // Edges use 0-based indices (converted in add_edge), so
            // valid range is [0, vocab_size).
            if a >= vocab_size {
                return Err(Error::AtomConnectionTokenError(
                    crate::errors::AtomConnectionTokenError::AtomIndexOutOfBounds {
                        index: a + 1,
                        num_atoms: vocab_size,
                    },
                ));
            }
            if b >= vocab_size {
                return Err(Error::AtomConnectionTokenError(
                    crate::errors::AtomConnectionTokenError::AtomIndexOutOfBounds {
                        index: b + 1,
                        num_atoms: vocab_size,
                    },
                ));
            }
        }

        edges.sort_unstable();
        edges.dedup();

        let edges: SymmetricCSR2D<CSR2D<Self::Idx, Self::Idx, Self::Idx>> = MolcularGraphEdgesBuilder::default()
            .expected_number_of_edges(Self::Idx::try_from(edges.len())?)
            .expected_shape(vocab_size)
            .edges(edges.clone().into_iter())
            .build()
            .unwrap_or_else(|_| panic!(
                "Failed to build UndiEdgesBuilder: vocabulary size: {vocab_size}, edges: {edges:?}, number of edges: {}",
                edges.len()
            ));

        Ok(MolecularGraph::from((vocab_size, edges)))
    }
}

impl PrefixFromStrWithContext for AtomConnectionLayer {}

/// Type alias for a generic undirected edges list builder.
pub type MolcularGraphEdgesBuilder<EdgeIterator, AtomIdx, EdgeIdx> =
    GenericUndirectedMonopartiteEdgesBuilder<
        EdgeIterator,
        UpperTriangularCSR2D<CSR2D<EdgeIdx, AtomIdx, AtomIdx>>,
        SymmetricCSR2D<CSR2D<EdgeIdx, AtomIdx, AtomIdx>>,
    >;
