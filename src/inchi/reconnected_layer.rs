//! Module for the reconnected layer of an InChI.

use crate::traits::prefix::Prefix;

#[derive(Debug, Clone, PartialEq, Eq)]
/// The reconnected layer of an InChI.
pub struct ReconnectedLayer;

impl Prefix for ReconnectedLayer {
    const PREFIX: char = 'r';
}
