//! Module defining the structure of the atom connection layer

use alloc::vec::Vec;

use geometric_traits::prelude::*;

use crate::traits::prefix::Prefix;
/// A molecular graph that is undirected.
pub type MolecularGraph<Idx> = GenericGraph<Idx, SymmetricCSR2D<CSR2D<Idx, Idx, Idx>>>;

/// The atom connection layer
pub(crate) type AtomConnectionLayer = Vec<MolecularGraph<usize>>;

impl Prefix for AtomConnectionLayer {
    const PREFIX: char = 'c';
}
