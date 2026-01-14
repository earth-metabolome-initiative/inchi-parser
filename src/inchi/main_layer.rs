use crate::traits::prefix::Prefix;

pub struct MainLayer {
    chemical_formula: ChemicalFormulaSubLayer,
    atom_connections: AtomConnectionsSubLayer,
    hydrogens: HydrogensSubLayer,
}

struct ChemicalFormulaSubLayer;
struct AtomConnectionsSubLayer;
struct HydrogensSubLayer;

impl Prefix for AtomConnectionsSubLayer {
    const PREFIX: char = 'c';
}

impl Prefix for HydrogensSubLayer {
    const PREFIX: char = 'h';
}
