//! Module for handling the charge layer.

use crate::traits::prefix::Prefix;

#[derive(Debug, Clone, PartialEq, Eq)]
/// Represents the charge layer of an InChI.
pub struct ChargeLayer {
    charged_sublayer: Option<ChargeSubLayer>,
    proton_sublayer: Option<ProtonSublayer>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ChargeSubLayer;
#[derive(Debug, Clone, PartialEq, Eq)]
struct ProtonSublayer;

impl Prefix for ChargeSubLayer {
    const PREFIX: char = 'q';
}

impl Prefix for ProtonSublayer {
    const PREFIX: char = 'p';
}
