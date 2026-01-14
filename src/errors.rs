//! Module for InChI-related errors.
use molecular_formulas::errors::Error as MolecularFormulaError;
use thiserror::Error;

/// Errors that can occur while parsing or handling InChIs.
#[derive(Error, Debug)]
pub enum Error {
    /// When the ""InChI="" prefix is missing.
    #[error("Missing InChI= prefix")]
    MissingInchiPrefix,
    /// When the version prefix is missing.
    #[error("Missing version prefix")]
    MissingVersionPrefix,
    /// When a layer prefix is missing a forward slash.
    #[error("Missing '/'")]
    MissingForwardSlash,
    /// Molecular formula errors
    #[error("Molecular formula error")]
    MolecularFormula(#[from] MolecularFormulaError),
    /// Is not Hill sorted error
    #[error("Molecular formula is not Hill sorted")]
    NotHillSorted,
    /// Wrong prefix for the layer/sublayer
    #[error("Wrong prefix for the layer")]
    WrongPrefix(String),
}
