//! Module for handling the charge layer.

use crate::traits::prefix::Prefix;

/// Represents the charge layer of an InChI.
pub struct ChargeLayer {
    charged_sublayer: Option<ChargeSubLayer>,
    proton_sublayer: Option<ProtonSublayer>,
}

struct ChargeSubLayer;
struct ProtonSublayer;

impl Prefix for ChargeSubLayer {
    const PREFIX: char = 'q';
}

impl Prefix for ProtonSublayer {
    const PREFIX: char = 'p';
}
