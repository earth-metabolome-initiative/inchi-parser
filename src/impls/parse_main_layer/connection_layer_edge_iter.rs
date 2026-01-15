//!

use std::str::Chars;

use crate::{
    errors::AtomConnectionEdgeError,
    impls::parse_main_layer::connection_layer_token_iter::{
        ConnectionLayerToken, ConnectionLayerTokenIter,
    },
};

pub(super) struct ConnectionLayerEdgeIter<I: Iterator<Item = char>> {
    token_iter: std::iter::Peekable<ConnectionLayerTokenIter<I>>,
    atom_index_stack: Vec<usize>,
}

impl<I: Iterator<Item = char>> From<ConnectionLayerTokenIter<I>> for ConnectionLayerEdgeIter<I> {
    fn from(token_iter: ConnectionLayerTokenIter<I>) -> Self {
        Self { token_iter: token_iter.peekable(), atom_index_stack: Vec::new() }
    }
}

impl<'a> From<&'a str> for ConnectionLayerEdgeIter<Chars<'a>> {
    fn from(value: &'a str) -> Self {
        let token_iter: ConnectionLayerTokenIter<_> = value.into();
        token_iter.into()
    }
}

impl<I: Iterator<Item = char>> Iterator for ConnectionLayerEdgeIter<I> {
    type Item = Result<(usize, usize), crate::errors::AtomConnectionEdgeError>;
    fn next(&mut self) -> Option<Self::Item> {
        let token = match self.token_iter.next()? {
            Ok(token) => token,
            Err(err) => return Some(Err(err.into())),
        };
        match token {
            ConnectionLayerToken::Atom(atom_index) => match self.token_iter.next() {
                Some(Err(err)) => return Some(Err(err.into())),
                None => return Some(Err(AtomConnectionEdgeError::OrphanAtomIndex(atom_index))),
                Some(Ok(ConnectionLayerToken::CloseRoundBracket)) => {
                    return Some(Err(AtomConnectionEdgeError::UnbalancedClosedParenthesis));
                }
                Some(Ok(ConnectionLayerToken::Atom(atom_index))) => {
                    unreachable!("An atom index cannot immediately follow an other atom index.")
                }
                Some(Ok(ConnectionLayerToken::OpenRoundBracket)) => match self.token_iter.next() {
                    Some(Ok(ConnectionLayerToken::Atom(atom_idex))) => {
                        todo!()
                    }
                    None => Some(Err(AtomConnectionEdgeError::UnbalancedOpenParenthesis)),
                    Some(Ok(unexpected_token)) => {
                        return Some(Err(
                            AtomConnectionEdgeError::UnexpectedTokenAfterParenthesis(
                                unexpected_token,
                            ),
                        ));
                    }
                    Some(Err(err)) => return Some(Err(err.into())),
                },
            },
        }
    }
}
