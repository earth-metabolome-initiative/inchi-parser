use core::str::FromStr;

use crate::{
    errors::Error,
    inchi::{InChI, MainLayer, proton_layer::ProtonSublayer},
    traits::{parse::ConsumeStr, prefix::Prefix},
    version::Version,
};

impl<V: Version> FromStr for InChI<V> {
    type Err = Error<u16>;

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
            return Err(Self::Err::MissingForwardSlash(
                "Missing chemical formula forward slash",
            ));
        };

        // If the next token is a lowercase 'p' then this means that we don't have
        // a main layer and we skip it
        if s.starts_with(ProtonSublayer::PREFIX) {
            return Err(Error::UnimplementedFeature(
                "Proton-only InChIs not yet supported",
            ));
        }

        // Parse the main layer (formula, connections, hydrogens).
        // Remaining layers after the hydrogen layer are validated but not parsed yet.
        let (main_layer, layer_remainder) = MainLayer::consume_str(s)?;

        // Validate that every remaining segment starts with a known layer prefix.
        if !layer_remainder.is_empty() {
            for segment in layer_remainder.split('/') {
                let Some(prefix) = segment.chars().next() else {
                    return Err(Error::UnrecognizedLayerPrefix('/'));
                };
                if !crate::constants::KNOWN_LAYER_PREFIXES.contains(&prefix) {
                    return Err(Error::UnrecognizedLayerPrefix(prefix));
                }
            }
        }

        Ok(InChI {
            main_layer,
            charge: None,
            stereochemistry: None,
            isotope: None,
            fixed_hydrogen: None,
            reconnected: None,
            _version: core::marker::PhantomData,
        })
    }
}
