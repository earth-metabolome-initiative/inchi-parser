//! Module for handling the charge layer.

use alloc::vec::Vec;

use crate::traits::prefix::Prefix;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ChargeSubLayer {
    pub(crate) charges: Vec<i16>,
}

impl Prefix for ChargeSubLayer {
    const PREFIX: char = 'q';
}
