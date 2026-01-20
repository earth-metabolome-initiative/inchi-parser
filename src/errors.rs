//! Module for InChI-related errors.
use crate::impls::parse_main_layer::connection_layer_base_token_iter::{
    ConnectionLayerSubToken, ConnectionLayerToken,
};
use core::num::NonZero;
use molecular_formulas::errors::Error as MolecularFormulaError;

/// Errors that can occur while parsing or handling InChIs.
#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum Error<Idx> {
    /// When the ""InChI="" prefix is missing.
    #[error("Missing InChI= prefix")]
    MissingInchiPrefix,
    /// When the version prefix is missing.
    #[error("Missing version prefix")]
    MissingVersionPrefix,
    /// When a layer prefix is missing a forward slash.
    #[error("Missing '/': {0}")]
    MissingForwardSlash(&'static str),
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
    AtomConnectionTokenError(#[from] AtomConnectionTokenError<Idx>),
    /// TODO! TEMPORARY ERROR TO REMOVE!
    #[error("Unimplemented feature: {0}")]
    UnimplementedFeature(&'static str),
}

/// Errors that can occur while tokenizing the atom connection layer.
#[derive(thiserror::Error, Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomConnectionTokenError<Idx> {
    /// Invalid token encountered
    #[error("Invalid character encountered: '{0}'")]
    InvalidCharacter(char),
    /// Invalid atom index encountered
    #[error("Atom index found larger than maximum size of the index type")]
    IndexOverflow,
    /// Zero atom index encountered
    #[error("Atom index cannot be zero")]
    IndexZero,
    /// The underlying iterator is intermittently empty
    #[error("The underlying iterator is empty")]
    UnderlyingIteratorEmpty,
    /// Two atom indices cannot be successive
    #[error("Two atom indices cannot be successive")]
    ConsecutiveAtomIndices,
    /// Illegal consecutive tokens
    #[error("Illegal consecutive tokens: '{previous}' followed by '{illegal}'")]
    IllegalConsecutiveSubTokens {
        /// The previous token
        previous: ConnectionLayerSubToken<Idx>,
        /// The illegal token
        illegal: ConnectionLayerSubToken<Idx>,
    },
    /// Unexpected end of input
    #[error("Unexpected end of input")]
    UnexpectedEndOfInput(ConnectionLayerSubToken<Idx>),
    /// Illegal closing bracket before an opening one
    #[error("There was a closing bracket before an open one.")]
    ClosingBracketBeforeOpeningBracket,
    /// If a comma appears before any edges
    #[error("Comma before any edge was added")]
    CommaBeforeAnyEdge,
    /// Illegal starting token
    #[error("Illegal starting token: '{0}'")]
    IllegalStartingToken(ConnectionLayerSubToken<Idx>),
    /// Self loop detected
    #[error("Self loop detected:  {0}")]
    SelfLoopDetected(Idx),
}
