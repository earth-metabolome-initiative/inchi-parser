use std::{fmt::Display, str::Chars};

#[derive(Debug)]
pub enum ConnectionLayerToken {
    OpenRoundBracket,
    CloseRoundBracket,
    Dash,
    Atom(usize),
    Comma,
}

impl Display for ConnectionLayerToken {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Dash => write!(f, "-"),
            Self::CloseRoundBracket => write!(f, ")"),
            Self::Comma => write!(f, ","),
            Self::OpenRoundBracket => write!(f, "("),
            Self::Atom(atom_index) => write!(f, "{atom_index}"),
        }
    }
}

/// Iterator over the `Token`s found in a provided string.
pub(super) struct ConnectionLayerTokenIter<I: Iterator<Item = char>> {
    /// The peekable chars iterator
    chars: std::iter::Peekable<I>,
}

impl<I> ConnectionLayerTokenIter<I>
where
    I: Iterator<Item = char>,
{
    fn parse_token(
        &mut self,
        current_char: char,
    ) -> Result<ConnectionLayerToken, crate::errors::AtomConnectionTokenError> {
        match current_char {
            '(' => Ok(ConnectionLayerToken::OpenRoundBracket),
            ')' => Ok(ConnectionLayerToken::CloseRoundBracket),
            '-' => Ok(ConnectionLayerToken::Dash),
            ',' => Ok(ConnectionLayerToken::Comma),
            c if c.is_ascii_digit() => {
                let mut number_str = String::new();
                number_str.push(c);
                while let Some(&next_char) = self.chars.peek() {
                    if next_char.is_ascii_digit() {
                        number_str.push(next_char);
                        self.chars.next();
                    } else {
                        break;
                    }
                }
                let number = number_str.parse::<usize>().map_err(|_| {
                    crate::errors::AtomConnectionTokenError::InvalidAtomIndex(number_str.clone())
                })?;
                Ok(ConnectionLayerToken::Atom(number))
            }
            _ => Err(crate::errors::AtomConnectionTokenError::InvalidToken(current_char)),
        }
    }
}

impl<I: Iterator<Item = char>> Iterator for ConnectionLayerTokenIter<I> {
    type Item = Result<ConnectionLayerToken, crate::errors::AtomConnectionTokenError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.chars.next().map(|current_char| self.parse_token(current_char))
    }
}

impl<I: Iterator<Item = char>> From<I> for ConnectionLayerTokenIter<I> {
    fn from(chars: I) -> Self {
        Self { chars: chars.peekable() }
    }
}

impl<'a> From<&'a str> for ConnectionLayerTokenIter<Chars<'a>> {
    fn from(value: &'a str) -> Self {
        value.chars().into()
    }
}
