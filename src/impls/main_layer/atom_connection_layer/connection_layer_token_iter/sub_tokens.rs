//! Submodule with the elemental tokens compositing a connection layer base
//! token.

use core::{fmt::Display, str::Chars};

use molecular_formulas::{BaselineDigit, NumberLike, try_fold_number};

use crate::{errors::AtomConnectionTokenError, traits::IndexLike};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
/// Enum representing the elemental tokens that compose a connection layer base
/// token.
pub enum ConnectionLayerSubToken<Idx> {
    /// Open parenthesis token: `(`.
    OpenParenthesis,
    /// Close parenthesis token: `)`.
    CloseParenthesis,
    /// Comma token: `,`.
    Comma,
    /// Dash token: `-`.
    Dash,
    /// Index token.
    Index(Idx),
}

impl<Idx> Display for ConnectionLayerSubToken<Idx>
where
    Idx: Display,
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ConnectionLayerSubToken::OpenParenthesis => write!(f, "("),
            ConnectionLayerSubToken::CloseParenthesis => write!(f, ")"),
            ConnectionLayerSubToken::Comma => write!(f, ","),
            ConnectionLayerSubToken::Dash => write!(f, "-"),
            ConnectionLayerSubToken::Index(idx) => write!(f, "{}", idx),
        }
    }
}

/// Iterator over the `Token`s found in a provided string.
pub(super) struct ConnectionLayerSubTokenIter<'a, Idx> {
    /// The peekable chars iterator
    chars: core::iter::Peekable<Chars<'a>>,
    /// Phantom data for the index type
    _phantom: core::marker::PhantomData<Idx>,
}

impl<'a, Idx> From<core::iter::Peekable<Chars<'a>>> for ConnectionLayerSubTokenIter<'a, Idx> {
    fn from(s: core::iter::Peekable<Chars<'a>>) -> Self {
        Self { chars: s, _phantom: core::marker::PhantomData }
    }
}

impl<'a, Idx: IndexLike> ConnectionLayerSubTokenIter<'a, Idx> {
    /// Returns whether the next character is a digit.
    pub fn peek_is_digit(&mut self) -> Option<bool> {
        Some(self.chars.peek()?.is_ascii_digit())
    }

    /// Returns whether the iterator is empty.
    pub fn is_empty(&mut self) -> bool {
        self.chars.peek().is_none()
    }
}

impl<Idx> Iterator for ConnectionLayerSubTokenIter<'_, Idx>
where
    Idx: IndexLike,
{
    type Item = Result<ConnectionLayerSubToken<Idx>, AtomConnectionTokenError<Idx>>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.peek_is_digit()? {
            let Some(maybe_number) = try_fold_number::<Idx, BaselineDigit, _>(&mut self.chars)
            else {
                unreachable!()
            };
            match maybe_number {
                Ok(idx) => return Some(Ok(ConnectionLayerSubToken::Index(idx))),
                Err(e) => return Some(Err(e.into())),
            }
        }
        let token = match self.chars.next()? {
            '(' => ConnectionLayerSubToken::OpenParenthesis,
            ')' => ConnectionLayerSubToken::CloseParenthesis,
            ',' => ConnectionLayerSubToken::Comma,
            '-' => ConnectionLayerSubToken::Dash,
            unexpected_char => {
                return Some(Err(AtomConnectionTokenError::InvalidCharacter(unexpected_char)));
            }
        };

        Some(match self.peek_is_digit() {
            Some(true) => {
                // If the next token is digit, we are done for this iteration.
                Ok(token)
            }
            Some(false) => {
                // If the next token is not a digit, we are in an error state.
                let next_token = match self.next()? {
                    Ok(next_token) => next_token,
                    Err(e) => return Some(Err(e)),
                };
                Err(AtomConnectionTokenError::IllegalConsecutiveSubTokens {
                    previous: token.into(),
                    illegal: next_token.into(),
                })
            }
            None => {
                // If the iterator is empty, we return an error as the current
                // token cannot be the last token.
                Err(AtomConnectionTokenError::UnexpectedEndOfInput(token))
            }
        })
    }
}
