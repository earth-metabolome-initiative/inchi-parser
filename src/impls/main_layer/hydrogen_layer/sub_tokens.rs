use core::{fmt::Display, str::Chars};

use molecular_formulas::{BaselineDigit, try_fold_number};

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
    /// The 'H' letter with an optional count (default 1).
    H(u8),
    /// Mobile group opener: `(H`, `(H2`, or `(H-`. Includes the consumed `(`, `H`,
    /// optional count, optional charge marker, and the mandatory `,` that follows.
    SharedHydrogens {
        /// Number of mobile hydrogens (1 if omitted).
        count: u8,
        /// `true` for `(H-,...)` charged mobile hydrogen (H⁻).
        charged: bool,
    },
    /// Close parenthesis token: `)`.
    CloseParenthesis,
    /// Asterisk token: `*`. Is always preceded by a digit.
    Asterisk(u8),
}

impl<Idx: Display> Display for HydogenLayerSubTokens<Idx> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            HydogenLayerSubTokens::SharedHydrogens { count, charged } => {
                write!(f, "(H")?;
                if *count != 1 {
                    write!(f, "{count}")?;
                }
                if *charged {
                    write!(f, "-")?;
                }
                Ok(())
            }
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

    /// Returns whether the next character is 'H'.
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
    u8: TryFrom<Idx, Error = core::num::TryFromIntError>,
{
    type Item = Result<HydogenLayerSubTokens<Idx>, HydrogenLayerTokenError<Idx>>;

    fn next(&mut self) -> Option<Self::Item> {
        // Handle digits first (atom index, range start, or repetition count)
        if let Some(true) = self.peek_is_digit() {
            let Some(maybe_number) = try_fold_number::<Idx, BaselineDigit, _>(&mut self.chars)
            else {
                unreachable!("digit confirmed by peek")
            };
            let n = match maybe_number {
                Ok(n) => n,
                Err(e) => return Some(Err(HydrogenLayerTokenError::NumericError(e))),
            };
            // Check what follows: '-' for range, '*' for repetition, else plain index
            if self.consume_dash() {
                // Range: N-M
                let m = match try_fold_number::<Idx, BaselineDigit, _>(&mut self.chars) {
                    Some(Ok(m)) => m,
                    Some(Err(e)) => return Some(Err(HydrogenLayerTokenError::NumericError(e))),
                    None => {
                        return Some(Err(HydrogenLayerTokenError::InvalidCharacter('-')));
                    }
                };
                return Some(Ok(HydogenLayerSubTokens::Range((n, m))));
            }
            if self.consume_asterisk() {
                return Some(match u8::try_from(n) {
                    Ok(count) => Ok(HydogenLayerSubTokens::Asterisk(count)),
                    Err(e) => Err(HydrogenLayerTokenError::TryFromIntError(e)),
                });
            }
            return Some(Ok(HydogenLayerSubTokens::Index(n)));
        }

        // Handle 'H' (fixed hydrogen marker, always outside parentheses at this point)
        if let Some(true) = self.peek_is_hydrogen() {
            self.chars.next(); // consume 'H'
            let count = match try_fold_number::<u8, BaselineDigit, _>(&mut self.chars) {
                Some(Ok(c)) => c,
                Some(Err(e)) => return Some(Err(HydrogenLayerTokenError::NumericError(e))),
                None => 1,
            };
            return Some(Ok(HydogenLayerSubTokens::H(count)));
        }

        // Handle other characters
        match self.chars.next()? {
            '(' => {
                // Expect 'H' immediately after '('
                match self.chars.next() {
                    Some('H') => {}
                    Some(c) => return Some(Err(HydrogenLayerTokenError::InvalidCharacter(c))),
                    None => {
                        return Some(Err(HydrogenLayerTokenError::InvalidCharacter('(')));
                    }
                }
                // Parse optional count (default 1)
                let count = match try_fold_number::<u8, BaselineDigit, _>(&mut self.chars) {
                    Some(Ok(c)) => c,
                    Some(Err(e)) => return Some(Err(HydrogenLayerTokenError::NumericError(e))),
                    None => 1,
                };
                // Check for '-' (charged H⁻ marker)
                let charged = if self.chars.peek() == Some(&'-') {
                    self.chars.next();
                    true
                } else {
                    false
                };
                // Expect ',' after the header
                match self.chars.next() {
                    Some(',') => {}
                    Some(c) => return Some(Err(HydrogenLayerTokenError::InvalidCharacter(c))),
                    None => {
                        return Some(Err(HydrogenLayerTokenError::InvalidCharacter('(')));
                    }
                }
                self.in_parenthesis = true;
                Some(Ok(HydogenLayerSubTokens::SharedHydrogens { count, charged }))
            }
            ')' => {
                if !self.in_parenthesis {
                    return Some(Err(HydrogenLayerTokenError::InvalidCharacter(')')));
                }
                self.in_parenthesis = false;
                Some(Ok(HydogenLayerSubTokens::CloseParenthesis))
            }
            ',' => Some(Ok(HydogenLayerSubTokens::Comma)),
            c => Some(Err(HydrogenLayerTokenError::InvalidCharacter(c))),
        }
    }
}
