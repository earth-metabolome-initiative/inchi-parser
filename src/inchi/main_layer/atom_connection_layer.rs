//! Module defining the structure of the atom connection layer

use alloc::{string::String, vec::Vec};

use geometric_traits::prelude::*;

use crate::traits::prefix::Prefix;
/// A molecular graph that is undirected.
pub type MolecularGraph<Idx> = GenericGraph<Idx, SymmetricCSR2D<CSR2D<Idx, Idx, Idx>>>;

/// A single component of the atom connection layer, holding both the parsed
/// graph and the original connection string for faithful round-trip display.
#[derive(Debug, Clone)]
pub(crate) struct AtomConnectionComponent {
    /// The parsed molecular graph.
    pub(crate) graph: MolecularGraph<u16>,
    /// The original connection body string (without prefix or repetition).
    pub(crate) connection_string: String,
}

impl PartialEq for AtomConnectionComponent {
    fn eq(&self, other: &Self) -> bool {
        self.graph == other.graph
    }
}

impl Eq for AtomConnectionComponent {}

/// The atom connection layer
pub(crate) type AtomConnectionLayer = Vec<AtomConnectionComponent>;

impl Prefix for AtomConnectionLayer {
    const PREFIX: char = 'c';
}
