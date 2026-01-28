//! Module for handling the charge layer.

use crate::traits::prefix::Prefix;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ChargeSubLayer;

impl Prefix for ChargeSubLayer {
    const PREFIX: char = 'q';
}
