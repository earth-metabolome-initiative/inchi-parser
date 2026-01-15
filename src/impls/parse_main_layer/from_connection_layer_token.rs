//! Parses a validated iterator of connection layer tokens into a vector of edges.

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

fn parse_atom_successor(
    atom_index: usize,
    edges: &mut Vec<(usize, usize)>,
    iter: &mut Peekable<ConnectionLayerTokenIter<std::str::Chars<'_>>>,
) -> Result<(), AtomConnectionTokenError> {
    if let Some(token) = iter.next().transpose()? {
        match token {
            ConnectionLayerToken::Dash => match iter.peek().copied().transpose()? {
                Some(ConnectionLayerToken::Atom(other_index)) => {
                    edges.push((atom_index, other_index));
                    Ok(())
                }
                Some(other_token) => Err(AtomConnectionTokenError::IllegalConsecutiveTokens(
                    ConnectionLayerToken::Dash,
                    other_token,
                )),
                None => {
                    Err(AtomConnectionTokenError::UnexpectedEndOfInput(ConnectionLayerToken::Dash))
                }
            },
            ConnectionLayerToken::OpenRoundBracket => match iter.peek().copied().transpose()? {
                Some(ConnectionLayerToken::Atom(other_index)) => {
                    edges.push((atom_index, other_index));
                    Ok(())
                }
                Some(other_token) => Err(AtomConnectionTokenError::IllegalConsecutiveTokens(
                    ConnectionLayerToken::Dash,
                    other_token,
                )),
                None => {
                    Err(AtomConnectionTokenError::UnexpectedEndOfInput(ConnectionLayerToken::Dash))
                }
            },
            other => unimplemented!("Unhandled token after Atom: {:?}", other),
        }
    } else {
        Ok(())
    }
}

impl FromConnectionLayer for Vec<(usize, usize)> {
    fn from_connection_layer_token(tokens: &str) -> Result<Self, AtomConnectionTokenError> {
        let mut edges = Vec::new();
        let iter: ConnectionLayerTokenIter<_> = tokens.into();
        let mut peekable = iter.peekable();

        while let Some(token) = peekable.next() {
            match token? {
                ConnectionLayerToken::Atom(atom_index) => {
                    parse_atom_successor(atom_index, &mut edges, &mut peekable)?
                }
                _ => unreachable!("Unhandled token type"),
            }
        }

        Ok(edges)
    }
}
