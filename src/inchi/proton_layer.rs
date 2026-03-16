//! Module defining the Proton sublayer of the InChI

use crate::traits::prefix::Prefix;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ProtonSublayer {
    pub(crate) proton_count: i16,
}

impl Prefix for ProtonSublayer {
    const PREFIX: char = 'p';
}
