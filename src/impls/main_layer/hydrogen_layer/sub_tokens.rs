use alloc::string::String;
use core::{fmt::Display, str::Chars};

use crate::{errors::HydrogenLayerTokenError, traits::IndexLike};
use molecular_formulas::{NumberLike, parsable::try_fold_number};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
/// Enum representing the elemental tokens that compose a hydrogen layer.
pub enum HydogenLayerSubTokens<Idx> {
    /// Index token.
    Index(Idx),
    /// A range of atoms (indices separated by a dash)
    Range((Idx, Idx)),
    /// Comma token: `,`.
    Comma,
    /// The 'H' letter
    H(u8),
    /// Open parenthesis token: `(`.
    SharedHydrogens(u8),
    /// Close parenthesis token: `)`.
    CloseParenthesis,
    /// Asterisk token : `*``. Is always preceded by a digit
    Asterisk(u8),
}

impl<Idx: Display> Display for HydogenLayerSubTokens<Idx> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            HydogenLayerSubTokens::SharedHydrogens(counter) => write!(f, "(H{counter}"),
            HydogenLayerSubTokens::CloseParenthesis => write!(f, ")"),
            HydogenLayerSubTokens::Comma => write!(f, ","),
            HydogenLayerSubTokens::Range((start_idx, end_idx)) => {
                write!(f, "{start_idx}-{end_idx}")
            }
            HydogenLayerSubTokens::Index(idx) => write!(f, "{idx}"),
            HydogenLayerSubTokens::H(count) => {
                if *count == 1 {
                    write!(f, "H")
                } else {
                    write!(f, "H{}", count)
                }
            }
            HydogenLayerSubTokens::Asterisk(counter) => write!(f, "{counter}*"),
        }
    }
}

pub(super) struct HydrogenLayerSubTokenIter<'a, Idx> {
    /// The peekable chars iterator
    chars: core::iter::Peekable<Chars<'a>>,
    /// If we are inside parenthesis
    in_parenthesis: bool,
    /// Phantom data for the index type
    _phantom: core::marker::PhantomData<Idx>,
}

impl<'a, Idx> From<&'a str> for HydrogenLayerSubTokenIter<'a, Idx> {
    fn from(s: &'a str) -> Self {
        Self { chars: s.chars().peekable(), _phantom: core::marker::PhantomData }
    }
}

impl<'a, Idx: IndexLike> HydrogenLayerSubTokenIter<'a, Idx> {
    /// Returns whether the next character is a digit.
    pub fn peek_is_digit(&mut self) -> Option<bool> {
        Some(self.chars.peek()?.is_ascii_digit())
    }

    fn consume_dash(&mut self) -> bool {
        if self.chars.peek() == Some(&'-') {
            self.chars.next();
            true
        } else {
            false
        }
    }

    fn consume_asterisk(&mut self) -> bool {
        if self.chars.peek() == Some(&'*') {
            self.chars.next();
            true
        } else {
            false
        }
    }
}

impl<Idx: IndexLike + NumberLike> Iterator for HydrogenLayerSubTokenIter<'_, Idx>
where
    u8: TryFrom<Idx>,
{
    type Item = Result<HydogenLayerSubTokens<Idx>, HydrogenLayerTokenError<Idx>>;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(maybe_index) = try_fold_number(&mut self.chars) {
            let begin_index = match maybe_index {
                Ok(begin_index) => begin_index,
                Err(e) => return Some(Err(e.into())),
            };

            if self.consume_dash() {
                if self.in_parenthesis {
                    return Some(Err(HydrogenLayerTokenError::InvalidCharacter('-')));
                }
                return Some(match try_fold_number(&mut self.chars) {
                    None => Err(HydrogenLayerTokenError::InvalidCharacter('-')),
                    Some(Err(e)) => Err(e.into()),
                    Some(Ok(end_index)) => {
                        Ok(HydogenLayerSubTokens::Range((begin_index, end_index)))
                    }
                });
            }

            if self.consume_asterisk() {
                if self.in_parenthesis {
                    return Some(Err(HydrogenLayerTokenError::InvalidCharacter('*')));
                }
                let Ok(index) = begin_index.try_into() else {
                    return Some(Err(HydrogenLayerTokenError::InvalidCharacter('*')));
                };
                return Some(Ok(HydogenLayerSubTokens::Asterisk(index)));
            }
        }

        let token = match self.chars.next()? {
            '(' => HydogenLayerSubTokens::OpenParenthesis,
            ')' => HydogenLayerSubTokens::CloseParenthesis,
            ',' => HydogenLayerSubTokens::Comma,
            'H' => {
                if self.peek_is_digit()? {
                    match self.consume_digit::<u8>() {
                        Ok(counter) => HydogenLayerSubTokens::H(counter),
                        Err(e) => return Some(Err(e)),
                    }
                } else {
                    HydogenLayerSubTokens::H(1)
                }
            }
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
