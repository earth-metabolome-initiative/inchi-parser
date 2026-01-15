//! Module for the main layer of an InChI.

use crate::traits::prefix::Prefix;
use geometric_traits::prelude::*;
use molecular_formulas::MolecularFormula;

/// The main layer of an InChI.
pub struct MainLayer {
    chemical_formula: MolecularFormula,
    atom_connections: Option<Vec<MolecularGraph>>,
    hydrogens: HydrogensSubLayer,
}

/// A molecular graph that is undirected.
pub type MolecularGraph = GenericGraph<usize, SymmetricCSR2D<CSR2D<usize, usize, usize>>>;

pub(crate) struct HydrogensSubLayer;

impl Prefix for Option<Vec<MolecularGraph>> {
    const PREFIX: char = 'c';
}

impl Prefix for HydrogensSubLayer {
    const PREFIX: char = 'h';
}
