use crate::inchi::InChI;
use crate::traits::parse::ParseLayer;
use crate::version::Version;
use molecular_formulas::MolecularFormula;
use std::str::FromStr;

impl<V: Version> FromStr for InChI<V> {
    type Err = crate::errors::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
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
            return Err(Self::Err::MissingForwardSlash);
        };

        // Then we parse the molecular formula layer
        // we strip everything until the next '/'
        let Some((mf_layer, rest)) = s.split_once('/') else {
            // if there is no '/' left, the InChI is invalid
            return Err(Self::Err::MissingForwardSlash);
        };

        // Then we parse the molecular formula layer
        let chemical_formula = MolecularFormula::parse(mf_layer, &mut ())?;

        // The atom connection layer
        todo!();

        // The hydrogen layer
        todo!();
    }
}
