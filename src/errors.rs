//! Module for InChI-related errors.
use crate::impls::parse_main_layer::connection_layer_base_token_iter::ConnectionLayerToken;
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
    /// TODO! TEMPORARY ERROR TO REMOVE!
    #[error("Unimplemented feature: {0}")]
    UnimplementedFeature(&'static str),
}

/// Errors that can occur while tokenizing the atom connection layer.
#[derive(thiserror::Error, Debug, Clone, Copy)]
pub enum AtomConnectionTokenError {
    /// Invalid token encountered
    #[error("Invalid token encountered: '{0}'")]
    InvalidToken(char),
    /// Invalid atom index encountered
    #[error("Atom index found larger than maximum size of {maximum_size}")]
    OverflowingAtomIndex {
        /// The maximum size allowed for atom indices
        maximum_size: usize,
    },
    /// The underlying iterator is intermittently empty
    #[error("The underlying iterator is empty")]
    UnderlyingIteratorEmpty,
    /// Two atom indices cannot be successive
    #[error("Two atom indices cannot be successive")]
    ConsecutiveAtomIndices,
    /// Illegal consecutive tokens
    #[error("Illegal consecutive tokens: '{0}' followed by '{1}'")]
    IllegalConsecutiveTokens(ConnectionLayerToken, ConnectionLayerToken),
    /// Unexpected end of input
    #[error("Unexpected end of input")]
    UnexpectedEndOfInput(ConnectionLayerToken),
}
