//! Parses a validated iterator of connection layer tokens into a vector of edges.
use core::num::NonZero;

use std::iter::Peekable;

use crate::{
    errors::AtomConnectionTokenError,
    impls::parse_main_layer::connection_layer_base_token_iter::{
        ConnectionLayerToken, ConnectionLayerTokenIter,
    },
};

pub(super) trait FromConnectionLayer: Sized {
    /// Parses a validated iterator of connection layer tokens into a vector of edges.
    ///
    /// # Arguments
    ///
    /// * `tokens` - A string slice representing the connection layer tokens.
    ///
    /// # Errors
    ///
    /// * Returns `AtomConnectionTokenError` if parsing fails.
    fn from_connection_layer_token(tokens: &str) -> Result<Self, AtomConnectionTokenError>;
}

/// Adds an edge between two atom indices to the edges vector.
fn add_edge(
    left_index: NonZero<usize>,
    right_index: NonZero<usize>,
    edges: &mut Vec<(usize, usize)>,
) -> Result<(), AtomConnectionTokenError> {
    if left_index == right_index {
        return Err(AtomConnectionTokenError::SelfLoopDetected(left_index));
    }

    let left = usize::from(left_index) - 1;
    let right = usize::from(right_index) - 1;
    edges.push(if left <= right { (left, right) } else { (right, left) });

    Ok(())
}

fn parse_atom_successor(
    atom_index: NonZero<usize>,
    edges: &mut Vec<(usize, usize)>,
    iter: &mut Peekable<ConnectionLayerTokenIter<'_>>,
) -> Result<(), AtomConnectionTokenError> {
    let Some(maybe_token) = iter.peek().copied().transpose()? else {
        return Ok(());
    };

    // If we encountered a closing bracket, we return to the previous level
    // and we do not consume the token here.
    if matches!(
        maybe_token.as_ref(),
        ConnectionLayerToken::CloseRoundBracket | ConnectionLayerToken::Comma
    ) {
        return Ok(());
    }

    // We consume the token here since we only peeked before.
    let Some(token) = iter.next().transpose()? else {
        return Ok(());
    };

    match token.as_ref() {
        ConnectionLayerToken::Atom(other_atom_index) => {
            add_edge(atom_index, *other_atom_index, edges)
        }
        ConnectionLayerToken::OpenRoundBracket => {
            let Some(next_token) = iter.next().transpose()? else {
                return Err(AtomConnectionTokenError::UnexpectedEndOfInput(token));
            };
            match next_token.as_ref() {
                ConnectionLayerToken::Atom(first_atom_in_brackets_index) => {
                    add_edge(atom_index, *first_atom_in_brackets_index, edges)?;
                    parse_atom_successor(*first_atom_in_brackets_index, edges, iter)?;
                    loop {
                        let Some(inner_token) = iter.next().transpose()? else {
                            todo!("AtomConnectionTokenError::UnexpectedEndOfInput");
                        };

                        return match inner_token.as_ref() {
                            ConnectionLayerToken::CloseRoundBracket => {
                                match iter.peek().copied().transpose()? {
                                    Some(ConnectionLayerToken::Atom(
                                        first_atom_out_of_brackets_index,
                                    )) => {
                                        add_edge(
                                            atom_index,
                                            first_atom_out_of_brackets_index,
                                            edges,
                                        )?;
                                        Ok(())
                                    }
                                    Some(other_token) => {
                                        Err(AtomConnectionTokenError::IllegalConsecutiveTokens(
                                            ConnectionLayerToken::CloseRoundBracket,
                                            other_token,
                                        ))
                                    }
                                    None => Err(AtomConnectionTokenError::UnexpectedEndOfInput(
                                        ConnectionLayerToken::CloseRoundBracket,
                                    )),
                                }
                            }
                            ConnectionLayerToken::Comma => match iter.next().transpose()? {
                                Some(ConnectionLayerToken::Atom(atom_index_after_comma)) => {
                                    add_edge(atom_index, atom_index_after_comma, edges)?;
                                    parse_atom_successor(atom_index_after_comma, edges, iter)?;
                                    continue;
                                }
                                Some(other_token) => {
                                    Err(AtomConnectionTokenError::IllegalConsecutiveTokens(
                                        ConnectionLayerToken::CloseRoundBracket,
                                        other_token,
                                    ))
                                }
                                None => Err(AtomConnectionTokenError::UnexpectedEndOfInput(
                                    ConnectionLayerToken::CloseRoundBracket,
                                )),
                            },
                            ConnectionLayerToken::Atom(_) => {
                                continue;
                            }
                            other_token => Err(AtomConnectionTokenError::IllegalConsecutiveTokens(
                                ConnectionLayerToken::OpenRoundBracket,
                                other_token,
                            )),
                        };
                    }
                }
                _ => Err(AtomConnectionTokenError::IllegalConsecutiveTokens(token, next_token)),
            }
        }
        other => unimplemented!("Unhandled token after Atom: {:?}", other),
    }
}

impl FromConnectionLayer for Vec<(usize, usize)> {
    fn from_connection_layer_token(tokens: &str) -> Result<Self, AtomConnectionTokenError> {
        let mut edges = Vec::new();
        let iter: ConnectionLayerTokenIter = tokens.into();
        let mut peekable = iter.peekable();

        while let Some(token) = peekable.next().transpose()? {
            match token.as_ref() {
                ConnectionLayerToken::Atom(atom_index) => {
                    parse_atom_successor(*atom_index, &mut edges, &mut peekable)?
                }
                ConnectionLayerToken::Comma
                | ConnectionLayerToken::OpenRoundBracket
                | ConnectionLayerToken::CloseRoundBracket => {
                    return Err(AtomConnectionTokenError::IllegalStartingToken(token));
                }
            }
        }

        Ok(edges)
    }
}
