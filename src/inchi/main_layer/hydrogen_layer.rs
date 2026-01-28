//!

use crate::traits::prefix::Prefix;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct HydrogensSubLayer;

impl Prefix for HydrogensSubLayer {
    const PREFIX: char = 'h';
}
