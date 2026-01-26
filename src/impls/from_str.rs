use core::str::FromStr;

use molecular_formulas::{InChIFormula, MolecularFormula};

use crate::{
    errors::Error,
    inchi::{
        InChI,
        main_layer::{AtomConnectionLayer, MolecularGraph},
    },
    traits::parse::ParseLayer,
    version::Version,
};

impl<V: Version> FromStr for InChI<V> {
    type Err = Error<usize>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // TODO: split the string at all '/'

        // First we remove the "InChI=" prefix
        let Some(s) = s.strip_prefix(crate::constants::INCHI_PREFIX) else {
            return Err(Self::Err::MissingInchiPrefix);
        };

        // Next the version prefix is removed
        let Some(s) = s.strip_prefix(V::VERSION_PREFIX) else {
            return Err(Self::Err::MissingVersionPrefix);
        };

        // Then we remove the first '/' to get the first layer
        let Some(s) = s.strip_prefix('/') else {
            return Err(Self::Err::MissingForwardSlash("Missing chemical formula forward slash"));
        };

        // Then we parse the molecular formula layer
        // we strip everything until the next '/'
        let Some((mf_layer, rest)) = s.split_once('/') else {
            // The hydrogen layer
            return Err(crate::errors::Error::UnimplementedFeature(
                "Hydrogen layer parsing not implemented",
            ));
        };

        // Then we parse the molecular formula layer
        let chemical_formula = InChIFormula::parse(mf_layer, ())?;

        let Some((next_layer, rest)) = rest.split_once('/') else {
            // The hydrogen layer
            return Err(crate::errors::Error::UnimplementedFeature(
                "Hydrogen layer parsing not implemented",
            ));
        };

        // The atom connection layer
        let undi_graphs = <AtomConnectionLayer>::parse(next_layer, &chemical_formula)?;

        // The hydrogen layer
        Err(crate::errors::Error::UnimplementedFeature("Hydrogen layer parsing not implemented"))
    }
}
