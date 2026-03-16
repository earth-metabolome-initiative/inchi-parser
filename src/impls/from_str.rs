use core::str::FromStr;

use crate::{
    errors::Error,
    inchi::{
        InChI, IsotopeLayer, MainLayer, charge_layer::ChargeSubLayer, proton_layer::ProtonSublayer,
    },
    traits::{
        parse::{ConsumeStr, PrefixFromStrWithContext},
        prefix::Prefix,
    },
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
            return Err(Self::Err::MissingForwardSlash("Missing chemical formula forward slash"));
        };

        // If the next token is a lowercase 'p' then this means that we don't have
        // a main layer and we skip it
        if s.starts_with(ProtonSublayer::PREFIX) {
            return Err(Error::UnimplementedFeature("Proton-only InChIs not yet supported"));
        }

        // Parse the main layer (formula, connections, hydrogens).
        // Remaining layers after the hydrogen layer are validated but not parsed yet.
        let (main_layer, mut layer_remainder) = MainLayer::consume_str(s)?;

        let charge =
            ChargeSubLayer::try_build_layer(&mut layer_remainder, main_layer.chemical_formula())?;

        let proton = ProtonSublayer::try_build_layer(&mut layer_remainder, ())?;

        // Skip stereo layers (b, t, m, s) that are not yet parsed.
        // These are known-good prefixes, so no further validation needed.
        while layer_remainder.starts_with('b')
            || layer_remainder.starts_with('t')
            || layer_remainder.starts_with('m')
            || layer_remainder.starts_with('s')
        {
            let (_, rest) = layer_remainder.split_once('/').unwrap_or(("", ""));
            layer_remainder = rest;
        }

        let isotope = IsotopeLayer::try_build_layer(
            &mut layer_remainder,
            (main_layer.chemical_formula(), None),
        )?;

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
            charge,
            proton,
            stereochemistry: None,
            isotope,
            fixed_hydrogen: None,
            reconnected: None,
            _version: core::marker::PhantomData,
        })
    }
}
