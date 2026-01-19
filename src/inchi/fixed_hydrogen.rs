//! Module for handling the fixed hydrogen layer.

use crate::traits::prefix::Prefix;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
/// Represents the fixed hydrogen layer of an InChI.
pub struct FixedHydrogenLayer;

impl Prefix for FixedHydrogenLayer {
    const PREFIX: char = 'f';
}
