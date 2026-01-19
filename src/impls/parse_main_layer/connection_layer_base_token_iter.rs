use core::num::NonZero;
use std::{fmt::Display, iter::Enumerate, str::Chars};
mod sub_tokens;
pub use sub_tokens::ConnectionLayerSubToken;
use sub_tokens::ConnectionLayerSubTokenIter;

use crate::errors::AtomConnectionTokenError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConnectionLayerToken {
    Branch(Box<ConnectionLayerToken>),
    Atom(NonZero<usize>),
    List(Vec<ConnectionLayerToken>),
}

impl Display for ConnectionLayerToken {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Branch(tokens) => {
                write!(f, "({})", tokens)
            }
            Self::Atom(atom_index) => write!(f, "{atom_index}"),
            Self::List(tokens) => {
                for (i, token) in tokens.iter().enumerate() {
                    if i > 0 {
                        write!(f, ",")?;
                    }
                    write!(f, "{token}")?;
                }
                Ok(())
            }
        }
    }
}

/// Iterator over the `Token`s found in a provided string.
pub(super) struct ConnectionLayerTokenIter<'a> {
    /// The peekable chars iterator
    tokens: ConnectionLayerSubTokenIter<'a>,
}

impl Iterator for ConnectionLayerTokenIter<'_> {
    type Item = Result<ConnectionLayerToken, crate::errors::AtomConnectionTokenError>;
    fn next(&mut self) -> Option<Self::Item> {
        Some(match self.next()? {
            '(' => {
                let Some(branch) = self.next().transpose()? else {
                    return Some(Err(AtomConnectionTokenError::UnclosedBranch(position)));
                };
            }
            ')' => Ok(ConnectionLayerToken::new_single_char(
                ConnectionLayerToken::CloseRoundBracket,
                position,
            )),
            '-' => return self.next(), // Skip hyphens
            ',' => Ok(ConnectionLayerToken::new_single_char(ConnectionLayerToken::Comma, position)),
            c if c.is_ascii_digit() => {
                let mut number_str = String::new();
                number_str.push(c);
                let mut current_position = position;
                while let Some(&(char_position, next_char)) = self.chars.peek() {
                    current_position = char_position;
                    if next_char.is_ascii_digit() {
                        number_str.push(next_char);
                        self.chars.next();
                    } else {
                        break;
                    }
                }
                let number = match number_str.parse::<NonZero<usize>>() {
                    Ok(n) => n,
                    Err(_) => {
                        return Some(Err(AtomConnectionTokenError::OverflowingAtomIndex {
                            maximum_size: usize::MAX,
                        }));
                    }
                };
                Ok(ConnectionLayerToken::new_atom(number, position, current_position))
            }
            _ => Err(AtomConnectionTokenError::InvalidToken(current_char)),
        })
    }
}

impl<'a> From<&'a str> for ConnectionLayerTokenIter<'a> {
    fn from(value: &'a str) -> Self {
        Self { chars: value.chars().enumerate().peekable() }
    }
}
