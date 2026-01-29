use alloc::string::String;
use core::{fmt::Display, str::Chars};

use crate::{errors::HydrogenLayerTokenError, traits::IndexLike};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
/// Enum representing the elemental tokens that compose a hydrogen layer.
pub enum HydogenLayerSubTokens<Idx> {
    /// Index token.
    Index(Idx),
    /// Dash token: `-`.
    Dash,
    /// Comma token: `,`.
    Comma,
    /// The 'H' letter
    H,
    /// Open parenthesis token: `(`.
    OpenParenthesis,
    /// Close parenthesis token: `)`.
    CloseParenthesis,
    /// Asterisk token : `*``
    Asterisk,
}

impl<Idx: Display> Display for HydogenLayerSubTokens<Idx> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            HydogenLayerSubTokens::OpenParenthesis => write!(f, "("),
            HydogenLayerSubTokens::CloseParenthesis => write!(f, ")"),
            HydogenLayerSubTokens::Comma => write!(f, ","),
            HydogenLayerSubTokens::Dash => write!(f, "-"),
            HydogenLayerSubTokens::Index(idx) => write!(f, "{}", idx),
            HydogenLayerSubTokens::H => write!(f, "H"),
            HydogenLayerSubTokens::Asterisk => write!(f, "*"),
        }
    }
}

pub(super) struct HydrogenLayerTokenIter<'a, Idx> {
    /// The peekable chars iterator
    chars: core::iter::Peekable<Chars<'a>>,
    /// Phantom data for the index type
    _phantom: core::marker::PhantomData<Idx>,
}

impl<'a, Idx> From<&'a str> for HydrogenLayerTokenIter<'a, Idx> {
    fn from(s: &'a str) -> Self {
        Self { chars: s.chars().peekable(), _phantom: core::marker::PhantomData }
    }
}

impl<'a, Idx: IndexLike> HydrogenLayerTokenIter<'a, Idx> {
    /// Returns whether the next character is a digit.
    pub fn peek_is_digit(&mut self) -> Option<bool> {
        Some(self.chars.peek()?.is_ascii_digit())
    }

    /// Consumes and returns a digit as a `Idx`
    ///
    /// # Errors
    ///
    /// * If the parsed digit overflows the `Idx` type.
    pub fn consume_digit(&mut self) -> Result<Idx, HydrogenLayerTokenError<Idx>> {
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
            return Err(HydrogenLayerTokenError::<Idx>::IndexOverflow);
        };

        if number == Idx::ZERO {
            return Err(HydrogenLayerTokenError::<Idx>::IndexOverflow);
        }

        Ok(number)
    }
}

impl<Idx: IndexLike> Iterator for HydrogenLayerTokenIter<'_, Idx> {
    type Item = Result<HydogenLayerSubTokens<Idx>, HydrogenLayerTokenError<Idx>>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.peek_is_digit()? {
            return match self.consume_digit() {
                Ok(idx) => Some(Ok(HydogenLayerSubTokens::Index(idx))),
                Err(e) => Some(Err(e)),
            };
        }
        let token = match self.chars.next()? {
            '(' => HydogenLayerSubTokens::OpenParenthesis,
            ')' => HydogenLayerSubTokens::CloseParenthesis,
            ',' => HydogenLayerSubTokens::Comma,
            '-' => HydogenLayerSubTokens::Dash,
            'H' => HydogenLayerSubTokens::H,
            '*' => HydogenLayerSubTokens::Asterisk,
            unexpected_char => {
                return Some(Err(HydrogenLayerTokenError::InvalidCharacter(unexpected_char)));
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
                Err(HydrogenLayerTokenError::IllegalConsecutiveSubTokens {
                    previous: token.into(),
                    illegal: next_token.into(),
                })
            }
            None => {
                // If the iterator is empty, we return an error as the current
                // token cannot be the last token.
                Err(HydrogenLayerTokenError::UnexpectedEndOfInput(token))
            }
        })
    }
}
