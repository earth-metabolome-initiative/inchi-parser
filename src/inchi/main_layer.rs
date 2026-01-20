//! Module for the main layer of an InChI.

use geometric_traits::prelude::*;
use molecular_formulas::MolecularFormula;

use crate::traits::prefix::Prefix;

/// The atom connection layer
pub type AtomConnectionLayer = Option<Vec<MolecularGraph<usize>>>;

#[derive(Debug, Clone, PartialEq, Eq)]
/// The main layer of an InChI.
pub struct MainLayer {
    chemical_formula: MolecularFormula,
    atom_connections: AtomConnectionLayer,
    hydrogens: HydrogensSubLayer,
}

/// A molecular graph that is undirected.
pub type MolecularGraph<Idx> = GenericGraph<Idx, SymmetricCSR2D<CSR2D<Idx, Idx, Idx>>>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct HydrogensSubLayer;

impl Prefix for AtomConnectionLayer {
    const PREFIX: char = 'c';
}

impl Prefix for HydrogensSubLayer {
    const PREFIX: char = 'h';
}
