use alloc::vec::Vec;
use core::fmt::Display;

use super::sub_tokens::HydrogenLayerSubTokenIter;

/// Iterator over the `Token`s found in a provided string.
pub(super) struct HydrogenLayerTokenIter<'a, Idx> {
    /// The peekable chars iterator
    tokens: HydrogenLayerSubTokenIter<'a, Idx>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HydrogenLayerToken<Idx> {
    /// An atom index
    Atom(Idx),
    Range,
    // Branch(Vec<Vec<HydrogenLayerToken<Idx>>>), // TODO: not sure how this should look like
}

// impl<Idx: IndexLike> Display for HydrogenLayerToken<Idx> {
//     fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
//         match self {
//             Self::Atom(atom_index) => write!(f, "{atom_index}"),
//             Self::Range => todo!(),
//         }
//     }
// }

impl<'a, Idx> From<&'a str> for HydrogenLayerTokenIter<'a, Idx> {
    fn from(value: &'a str) -> Self {
        Self { tokens: HydrogenLayerSubTokenIter::from(value) }
    }
}
