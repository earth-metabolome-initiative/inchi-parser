//! Module for InChI-related errors.
use crate::impls::parse_main_layer::connection_layer_token_iter::ConnectionLayerToken;
use geometric_traits::errors;
use molecular_formulas::errors::Error as MolecularFormulaError;

/// Errors that can occur while parsing or handling InChIs.
#[derive(thiserror::Error, Debug)]
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
    WrongPrefix,
    /// If the molecular formula contains a mixture but the atom connections layer does not
    #[error("Molecular formula contains {0} mixtures but atom connection has {1}.")]
    FormulaAndConnectionLayerMixtureMismatch(usize, usize),
    /// Errors while tokenizing the atom connection layer
    #[error("Atom connection tokenization error: {0}")]
    AtomConnectionTokenError(#[from] AtomConnectionTokenError),
}

/// Errors that can occur while tokenizing the atom connection layer.
#[derive(thiserror::Error, Debug)]
pub enum AtomConnectionTokenError {
    /// Invalid token encountered
    #[error("Invalid token encountered: '{0}'")]
    InvalidToken(char),
    /// Invalid atom index encountered
    #[error("Invalid atom index encountered: '{0}'")]
    InvalidAtomIndex(String),
}

#[derive(thiserror::Error, Debug)]
pub enum AtomConnectionEdgeError {
    /// A token error in the atom connection layer
    #[error("Token error")]
    AtomConnectionTokenError(#[from] AtomConnectionTokenError),
    /// Orphan atom index : The atom is present in the atom connection layer but has no edges to attach to
    #[error("The atom {0} is present in the atom connection layer but has no edges.")]
    OrphanAtomIndex(usize),
    /// Unbalanced closed parenthesis
    #[error("Unbalanced closed parenthesis")]
    UnbalancedClosedParenthesis,
    /// Unbalanced open parenthesis
    #[error("Unbalanced open parenthesis")]
    UnbalancedOpenParenthesis,
    /// Unexpected token after open parenthesis
    #[error("An unexpected toekn after parenthesis : '{0}'")]
    UnexpectedTokenAfterParenthesis(ConnectionLayerToken),
}
