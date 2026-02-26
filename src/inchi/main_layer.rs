//! Module for the main layer of an InChI.

mod atom_connection_layer;
pub(crate) use atom_connection_layer::{AtomConnectionLayer, MolecularGraph};
mod hydrogen_layer;
use core::str::FromStr;

pub(crate) use hydrogen_layer::{HydrogenComponent, HydrogensSubLayer, MobileHydrogenGroup};
use molecular_formulas::InChIFormula;

use crate::traits::parse::{ConsumeStr, PrefixFromStrWithContext};

#[derive(Debug, Clone, PartialEq, Eq)]
/// The main layer of an InChI.
pub struct MainLayer {
    chemical_formula: InChIFormula,
    atom_connections: Option<AtomConnectionLayer<u16>>,
    hydrogens: Option<HydrogensSubLayer>,
}

impl ConsumeStr for MainLayer {
    type Idx = u16;
    fn consume_str(input: &str) -> Result<(Self, &str), crate::errors::Error<Self::Idx>> {
        // Then we parse the molecular formula layer
        // we strip everything until the next '/'
        let (chemical_formula_layer, mut layer_remainder) =
            input.split_once('/').unwrap_or((input, ""));

        // Then we parse the molecular formula layer
        let chemical_formula = InChIFormula::from_str(chemical_formula_layer)?;

        let atom_connection_layer =
            AtomConnectionLayer::try_build_layer(&mut layer_remainder, &chemical_formula)?;

        let hydrogen_layer =
            HydrogensSubLayer::try_build_layer(&mut layer_remainder, &chemical_formula)?;

        Ok((
            MainLayer {
                chemical_formula: chemical_formula,
                atom_connections: atom_connection_layer,
                hydrogens: hydrogen_layer,
            },
            layer_remainder,
        ))
    }
}
