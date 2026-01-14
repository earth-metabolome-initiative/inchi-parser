//! Module for the reconnected layer of an InChI.

use crate::traits::prefix::Prefix;

/// The reconnected layer of an InChI.
pub struct ReconnectedLayer;

impl Prefix for ReconnectedLayer {
    const PREFIX: char = 'r';
}
