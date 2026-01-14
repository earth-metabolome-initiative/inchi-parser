//! Module for the isotope layer of an InChI.

use crate::traits::prefix::Prefix;

/// The isotope layer of an InChI.
pub struct IsotopeLayer;

impl Prefix for IsotopeLayer {
    const PREFIX: char = 'i';
}
