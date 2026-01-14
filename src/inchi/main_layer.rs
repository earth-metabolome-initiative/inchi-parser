//! Module for the main layer of an InChI.

use crate::traits::prefix::Prefix;
use molecular_formulas::MolecularFormula;

/// The main layer of an InChI.
pub struct MainLayer {
    chemical_formula: MolecularFormula,
    atom_connections: AtomConnectionsSubLayer,
    hydrogens: HydrogensSubLayer,
}

pub(crate) struct AtomConnectionsSubLayer;
pub(crate) struct HydrogensSubLayer;

impl Prefix for AtomConnectionsSubLayer {
    const PREFIX: char = 'c';
}

impl Prefix for HydrogensSubLayer {
    const PREFIX: char = 'h';
}
