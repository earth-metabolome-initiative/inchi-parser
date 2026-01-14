//! Module for handling the fixed hydrogen layer.

use crate::traits::prefix::Prefix;

/// Represents the fixed hydrogen layer of an InChI.
pub struct FixedHydrogenLayer;

impl Prefix for FixedHydrogenLayer {
    const PREFIX: char = 'f';
}
