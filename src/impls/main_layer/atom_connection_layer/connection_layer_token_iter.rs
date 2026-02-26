use core::{fmt::Display, str::Chars};
mod sub_tokens;
use alloc::vec::Vec;

pub use sub_tokens::ConnectionLayerSubToken;
use sub_tokens::ConnectionLayerSubTokenIter;

use crate::{errors::AtomConnectionTokenError, traits::IndexLike};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConnectionLayerToken<Idx> {
    Branch(Vec<Vec<ConnectionLayerToken<Idx>>>),
    Atom(Idx),
}

impl<Idx: IndexLike> Display for ConnectionLayerToken<Idx> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Branch(tokens) => {
                let mut first = true;
                write!(f, "(")?;
                for subtokens in tokens {
                    if !first {
                        write!(f, ",")?;
                    }
                    let mut last_was_atom = false;
                    for token in subtokens {
                        if last_was_atom && matches!(token, ConnectionLayerToken::Atom(_)) {
                            write!(f, "-")?;
                        }
                        last_was_atom = matches!(token, ConnectionLayerToken::Atom(_));
                        write!(f, "{token}")?;
                    }
                    first = false;
                }
                write!(f, ")")
            }
            Self::Atom(atom_index) => write!(f, "{atom_index}"),
        }
    }
}

/// Iterator over the `Token`s found in a provided string.
pub(super) struct ConnectionLayerTokenIter<'a, Idx> {
    /// The peekable chars iterator
    tokens: ConnectionLayerSubTokenIter<'a, Idx>,
}

impl<Idx: IndexLike> Iterator for ConnectionLayerTokenIter<'_, Idx> {
    type Item = Result<ConnectionLayerToken<Idx>, AtomConnectionTokenError<Idx>>;
    fn next(&mut self) -> Option<Self::Item> {
        let token = match self.tokens.next()? {
            Ok(sub_token) => sub_token,
            Err(e) => return Some(Err(e)),
        };
        Some(Ok(match token {
            ConnectionLayerSubToken::Index(atom_index) => ConnectionLayerToken::Atom(atom_index),
            // The dash can be ignored as we have already validated the illegal
            // states associated with it in the sub-token iterator.
            ConnectionLayerSubToken::Dash => return self.next(),
            ConnectionLayerSubToken::OpenParenthesis => {
                // We start constructing a branch token, which
                // might contain one or more sub-tokens.
                let mut branch_tokens = Vec::new();
                let mut sub_tokens = Vec::new();
                loop {
                    let sub_token = match self.next()? {
                        Ok(sub_token) => sub_token,
                        Err(e) => match e {
                            AtomConnectionTokenError::ClosingBracketBeforeOpeningBracket => {
                                break;
                            }
                            AtomConnectionTokenError::CommaBeforeAnyEdge => {
                                if sub_tokens.is_empty() {
                                    return Some(Err(e));
                                }
                                branch_tokens.push(sub_tokens);
                                sub_tokens = Vec::new();
                                continue;
                            }
                            _ => return Some(Err(e)),
                        },
                    };
                    sub_tokens.push(sub_token);
                }
                ConnectionLayerToken::Branch(branch_tokens)
            }
            ConnectionLayerSubToken::CloseParenthesis => {
                return Some(Err(AtomConnectionTokenError::ClosingBracketBeforeOpeningBracket));
            }
            ConnectionLayerSubToken::Comma => {
                return Some(Err(AtomConnectionTokenError::CommaBeforeAnyEdge));
            }
        }))
    }
}

impl<'a, Idx> From<core::iter::Peekable<Chars<'a>>> for ConnectionLayerTokenIter<'a, Idx> {
    fn from(value: core::iter::Peekable<Chars<'a>>) -> Self {
        Self { tokens: ConnectionLayerSubTokenIter::from(value) }
    }
}
