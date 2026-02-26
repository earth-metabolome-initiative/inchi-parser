//! Parses a validated iterator of connection layer tokens into a vector of
//! edges.

use core::str::Chars;

use alloc::vec::Vec;
use molecular_formulas::NumberLike;

use crate::{
    errors::AtomConnectionTokenError,
    impls::main_layer::atom_connection_layer::connection_layer_token_iter::{
        ConnectionLayerToken, ConnectionLayerTokenIter,
    },
    traits::IndexLike,
};

pub(super) trait FromConnectionLayer: Sized {
    /// Index type used in the edges.
    type AtomIndex: IndexLike;

    /// Parses a validated iterator of connection layer tokens into a vector of
    /// edges.
    ///
    /// # Arguments
    ///
    /// * `tokens` - A string slice representing the connection layer tokens.
    ///
    /// # Errors
    ///
    /// * Returns `AtomConnectionTokenError` if parsing fails.
    fn from_connection_layer_token(
        tokens: core::iter::Peekable<Chars<'_>>,
    ) -> Result<Self, AtomConnectionTokenError<Self::AtomIndex>>;
}

/// Adds an edge between two atom indices to the edges vector.
fn add_edge<Idx: IndexLike>(
    left_index: Idx,
    right_index: Idx,
    edges: &mut Vec<(Idx, Idx)>,
) -> Result<(), AtomConnectionTokenError<Idx>> {
    if left_index == right_index {
        return Err(AtomConnectionTokenError::SelfLoopDetected(left_index));
    }

    let left = left_index - Idx::ONE;
    let right = right_index - Idx::ONE;
    edges.push(if left <= right { (left, right) } else { (right, left) });

    Ok(())
}

fn parse_token<Idx: IndexLike>(
    last_atom: Option<Idx>,
    token: ConnectionLayerToken<Idx>,
    edges: &mut Vec<(Idx, Idx)>,
) -> Result<Idx, AtomConnectionTokenError<Idx>> {
    Ok(match token {
        ConnectionLayerToken::Atom(current_atom) => {
            if let Some(last) = last_atom {
                add_edge(last, current_atom, edges)?;
            }
            current_atom
        }
        ConnectionLayerToken::Branch(branch) => {
            let Some(last_atom) = last_atom else {
                unreachable!("Branch cannot be the first token")
            };
            for sub_branch in branch {
                let mut branch_last_atom = last_atom;
                for sub_token in sub_branch {
                    branch_last_atom = parse_token(Some(branch_last_atom), sub_token, edges)?;
                }
            }
            last_atom
        }
    })
}

impl<Idx: IndexLike> FromConnectionLayer for Vec<(Idx, Idx)> {
    type AtomIndex = Idx;

    fn from_connection_layer_token(
        tokens: core::iter::Peekable<Chars<'_>>,
    ) -> Result<Self, AtomConnectionTokenError<Self::AtomIndex>> {
        let mut edges = Vec::new();
        let iter: ConnectionLayerTokenIter<Idx> = tokens.into();
        let mut last_atom: Option<Idx> = None;
        for maybe_token in iter {
            let token = maybe_token?;
            last_atom = Some(parse_token(last_atom, token, &mut edges)?);
        }

        Ok(edges)
    }
}
