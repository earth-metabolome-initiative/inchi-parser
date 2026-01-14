use crate::traits::prefix::Prefix;
use crate::{inchi::main_layer::AtomConnectionsSubLayer, traits::parse::ParseLayer};
use molecular_formulas::MolecularFormula;
use std::str::FromStr;

impl ParseLayer for MolecularFormula {
    type Context = ();
    type Error = crate::errors::Error;
    fn parse(input: &str, _context: &mut Self::Context) -> Result<Self, Self::Error> {
        // parse the molecular formula
        let molecular_formula = MolecularFormula::from_str(input)?;
        let is_hill_sorted = molecular_formula.is_hill_sorted()?;
        if !is_hill_sorted {
            return Err(Self::Error::NotHillSorted);
        }
        Ok(molecular_formula)
    }
}

impl ParseLayer for AtomConnectionsSubLayer {
    type Context = MolecularFormula;
    type Error = crate::errors::Error;
    fn parse(input: &str, context: &mut Self::Context) -> Result<Self, Self::Error> {
        let Some(s) = input.strip_prefix(Self::PREFIX) else {
            return Err(Self::Error::WrongPrefix);
        };
    }
}
