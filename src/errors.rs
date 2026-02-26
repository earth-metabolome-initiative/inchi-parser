//! Module for InChI-related errors.

use molecular_formulas::errors::{NumericError, ParserError};

use crate::impls::main_layer::{
    atom_connection_layer::connection_layer_token_iter::ConnectionLayerSubToken,
    hydrogen_layer::sub_tokens::HydogenLayerSubTokens,
};

/// Errors that can occur while parsing or handling InChIs.
#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum Error<Idx> {
    /// When the ""InChI="" prefix is missing.
    #[error("Missing 'InChI=' prefix")]
    MissingInchiPrefix,
    /// When the version prefix is missing.
    #[error("Missing version prefix")]
    MissingVersionPrefix,
    /// When a layer prefix is missing a forward slash.
    #[error("Missing '/': {0}")]
    MissingForwardSlash(&'static str),
    /// Molecular formula parse errors
    #[error("Molecular formula parsing error")]
    MolecularFormulaParserError(#[from] ParserError),
    /// Wrong prefix for the layer/sublayer
    #[error("Wrong prefix for the layer")]
    WrongPrefix,
    /// If the molecular formula contains a mixture but the atom connections
    /// layer does not
    #[error("Molecular formula contains {0} mixtures but atom connection layer does not match.")]
    FormulaAndConnectionLayerMixtureMismatch(usize),
    /// Errors while tokenizing the atom connection layer
    #[error("Atom connection tokenization error: {0}")]
    AtomConnectionTokenError(#[from] AtomConnectionTokenError<Idx>),
    /// Errors when converting a numberic value to an other numeric format
    #[error("Numberic error: {0}")]
    TryFromIntError(#[from] core::num::TryFromIntError),
    /// Errors while tokenizing the hydrogen layer
    #[error("Hydrogen layer tokenization error: {0}")]
    HydrogenLayerTokenError(#[from] HydrogenLayerTokenError<Idx>),
    /// TODO! TEMPORARY ERROR TO REMOVE!
    #[error("Unimplemented feature: {0}")]
    UnimplementedFeature(&'static str),
}

/// Errors that can occur while tokenizing the atom connection layer.
#[derive(thiserror::Error, Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomConnectionTokenError<Idx> {
    /// Numeric error
    #[error("Atom index cannot be zero")]
    NumericError(#[from] NumericError),
    /// Invalid token encountered
    #[error("Invalid character encountered: '{0}'")]
    InvalidCharacter(char),
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

#[derive(thiserror::Error, Debug, Clone, Copy, PartialEq, Eq)]
/// Errors that can occur while tokenizing the hydrogen layer.
pub enum HydrogenLayerTokenError<Idx> {
    /// Numeric error
    #[error("Atom index cannot be zero")]
    NumericError(#[from] NumericError),
    /// Invalid token encountered
    #[error("Invalid character encountered: '{0}'")]
    InvalidCharacter(char),
    /// Illegal consecutive tokens
    #[error("Illegal consecutive tokens: '{previous}' followed by '{illegal}'")]
    IllegalConsecutiveSubTokens {
        /// The previous token
        previous: HydogenLayerSubTokens<Idx>,
        /// The illegal token
        illegal: HydogenLayerSubTokens<Idx>,
    },
    /// Unexpected end of input
    #[error("Unexpected end of input")]
    UnexpectedEndOfInput(HydogenLayerSubTokens<Idx>),
    /// Integer conversion error (e.g. hydrogen count too large for u8)
    #[error("Integer conversion error: {0}")]
    TryFromIntError(core::num::TryFromIntError),
}
