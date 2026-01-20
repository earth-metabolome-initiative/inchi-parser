//! Submodule with the elemental tokens compositing a connection layer base
//! token.

use core::{fmt::Display, str::Chars};

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
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
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
    chars: std::iter::Peekable<Chars<'a>>,
    /// Phantom data for the index type
    _phantom: std::marker::PhantomData<Idx>,
}

impl<'a, Idx> From<&'a str> for ConnectionLayerSubTokenIter<'a, Idx> {
    fn from(s: &'a str) -> Self {
        Self { chars: s.chars().peekable(), _phantom: std::marker::PhantomData }
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

    /// Consumes and returns a digit as a `Idx`
    ///
    /// # Errors
    ///
    /// * If the parsed digit overflows the `Idx` type.
    pub fn consume_digit(&mut self) -> Result<Idx, AtomConnectionTokenError<Idx>> {
        let mut number_str = String::new();
        while let Some(&next_char) = self.chars.peek() {
            if next_char.is_ascii_digit() {
                number_str.push(next_char);
                self.chars.next();
            } else {
                break;
            }
        }
        let Ok(number) = number_str.parse::<Idx>() else {
            return Err(AtomConnectionTokenError::IndexOverflow);
        };

        if number == Idx::ZERO {
            return Err(AtomConnectionTokenError::IndexZero);
        }

        Ok(number)
    }
}

impl<Idx: IndexLike> Iterator for ConnectionLayerSubTokenIter<'_, Idx> {
    type Item = Result<ConnectionLayerSubToken<Idx>, AtomConnectionTokenError<Idx>>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.peek_is_digit()? {
            return match self.consume_digit() {
                Ok(idx) => Some(Ok(ConnectionLayerSubToken::Index(idx))),
                Err(e) => Some(Err(e)),
            };
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
