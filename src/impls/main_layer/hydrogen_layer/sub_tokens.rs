use core::{fmt::Display, str::Chars};

use molecular_formulas::{NumberLike, parsable::try_fold_number};

use crate::{errors::HydrogenLayerTokenError, traits::IndexLike};

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
    /// Asterisk token : `*`. Is always preceded by a digit
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
        Self {
            chars: s.chars().peekable(),
            _phantom: core::marker::PhantomData,
            in_parenthesis: false,
        }
    }
}

impl<'a, Idx: IndexLike> HydrogenLayerSubTokenIter<'a, Idx> {
    /// Returns whether the next character is a digit.
    pub fn peek_is_digit(&mut self) -> Option<bool> {
        Some(self.chars.peek()?.is_ascii_digit())
    }

    pub fn peek_is_hydrogen(&mut self) -> Option<bool> {
        Some(self.chars.peek()? == &'H')
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

impl<Idx: IndexLike> Iterator for HydrogenLayerSubTokenIter<'_, Idx>
where
    u8: TryFrom<Idx>,
{
    type Item = Result<HydogenLayerSubTokens<Idx>, HydrogenLayerTokenError<Idx>>;
    fn next(&mut self) -> Option<Self::Item> {
        todo!()
    }
}
