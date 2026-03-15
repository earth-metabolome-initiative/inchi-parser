use alloc::vec::Vec;
use core::str::FromStr;

use crate::{
    errors::Error,
    inchi::{
        InChI, IsotopeLayer, MainLayer,
        charge_layer::ChargeSubLayer,
        isotope_layer::IsotopeComponent,
        proton_layer::ProtonSublayer,
        stereochemistry_layer::{
            AlleneSublayer, DoubleBondSublayer, StereoChemistryInformationSublayer,
            StereochemistryLayer, TetrahedralSublayer,
        },
    },
    traits::{
        parse::{ConsumeStr, PrefixFromStrWithContext},
        prefix::Prefix,
    },
    version::Version,
};

/// Parses a proton-only InChI (no chemical formula).
///
/// After the `/` separator, the input starts with `p` and may optionally
/// include an isotope layer in the form `/i/hXY`.
fn parse_proton_only<V: Version>(mut s: &str) -> Result<InChI<V>, Error<u16>> {
    let proton = ProtonSublayer::try_build_layer(&mut s, ())?;

    // Parse optional isotope layer: only `/i/hXY` form is valid here
    // (empty atom body + hydrogen isotope sublayer).
    let isotope = if s.starts_with(IsotopeLayer::PREFIX) {
        // Consume the "i" segment
        let (_, mut remainder) = s.split_once('/').unwrap_or((s, ""));

        // Check for hydrogen isotope sublayer (/hD, /hT, /hH)
        let hydrogens = if remainder.starts_with('h')
            && remainder.as_bytes().get(1).is_some_and(|&b| b == b'D' || b == b'T' || b == b'H')
        {
            let (h_seg, rest) = remainder.split_once('/').unwrap_or((remainder, ""));
            remainder = rest;
            super::isotope_layer::parse_h_isotope_segment(h_seg)?
        } else {
            Vec::new()
        };

        s = remainder;
        Some(IsotopeLayer {
            components: alloc::vec![IsotopeComponent { atoms: Vec::new(), hydrogens }],
        })
    } else {
        None
    };

    if !s.is_empty() {
        return Err(Error::UnrecognizedLayerPrefix(s.chars().next().unwrap_or('/')));
    }

    Ok(InChI {
        main_layer: None,
        charge: None,
        proton,
        stereochemistry: None,
        isotope,
        fixed_hydrogen: None,
        reconnected: None,
        _version: core::marker::PhantomData,
    })
}

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

        // Proton-only InChIs have no chemical formula — just /p and optionally /i
        if s.starts_with(ProtonSublayer::PREFIX) {
            return parse_proton_only(s);
        }

        // Parse the main layer (formula, connections, hydrogens).
        // Remaining layers after the hydrogen layer are validated but not parsed yet.
        let (main_layer, mut layer_remainder) = MainLayer::consume_str(s)?;

        let charge =
            ChargeSubLayer::try_build_layer(&mut layer_remainder, main_layer.chemical_formula())?;

        let proton = ProtonSublayer::try_build_layer(&mut layer_remainder, ())?;

        let double_bond = DoubleBondSublayer::try_build_layer(
            &mut layer_remainder,
            main_layer.chemical_formula(),
        )?;
        let tetrahedral = TetrahedralSublayer::try_build_layer(
            &mut layer_remainder,
            main_layer.chemical_formula(),
        )?;
        let allene = AlleneSublayer::try_build_layer(&mut layer_remainder, ())?;
        let stereo_info =
            StereoChemistryInformationSublayer::try_build_layer(&mut layer_remainder, ())?;

        let stereochemistry = if double_bond.is_some()
            || tetrahedral.is_some()
            || allene.is_some()
            || stereo_info.is_some()
        {
            Some(StereochemistryLayer { double_bond, tetrahedral, allene, stereo_info })
        } else {
            None
        };

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
            main_layer: Some(main_layer),
            charge,
            proton,
            stereochemistry,
            isotope,
            fixed_hydrogen: None,
            reconnected: None,
            _version: core::marker::PhantomData,
        })
    }
}
