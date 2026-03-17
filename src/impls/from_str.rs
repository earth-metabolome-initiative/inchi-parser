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

/// Four optional stereo sublayers parsed as a group.
type StereoSublayers = (
    Option<DoubleBondSublayer>,
    Option<TetrahedralSublayer>,
    Option<AlleneSublayer>,
    Option<StereoChemistryInformationSublayer>,
);

/// Try to parse one set of stereo sublayers (`/b`, `/t`, `/m`, `/s`) from
/// the front of `input`, returning each as an `Option`.
fn parse_stereo_sublayers(
    input: &mut &str,
    formula: &molecular_formulas::InChIFormula,
) -> Result<StereoSublayers, Error<u16>> {
    let db = DoubleBondSublayer::try_build_layer(input, formula)?;
    let tet = TetrahedralSublayer::try_build_layer(input, formula)?;
    let allene = AlleneSublayer::try_build_layer(input, ())?;
    let info = StereoChemistryInformationSublayer::try_build_layer(input, ())?;
    Ok((db, tet, allene, info))
}

/// Build a [`StereochemistryLayer`] from individual `Option` sublayers,
/// returning `None` when all sublayers are absent.
fn build_stereo_layer(
    double_bond: Option<DoubleBondSublayer>,
    tetrahedral: Option<TetrahedralSublayer>,
    allene: Option<AlleneSublayer>,
    stereo_info: Option<StereoChemistryInformationSublayer>,
) -> Option<StereochemistryLayer> {
    if double_bond.is_some() || tetrahedral.is_some() || allene.is_some() || stereo_info.is_some() {
        Some(StereochemistryLayer { double_bond, tetrahedral, allene, stereo_info })
    } else {
        None
    }
}

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
        isotope_stereochemistry: None,
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

        let (mut db, mut tet, mut allene, mut s_info) =
            parse_stereo_sublayers(&mut layer_remainder, main_layer.chemical_formula())?;

        let isotope = IsotopeLayer::try_build_layer(
            &mut layer_remainder,
            (main_layer.chemical_formula(), None),
        )?;

        // Parse stereo layers after /i.  Two scenarios produce these:
        //
        // 1. Non-standard ordering — /i appears before /b, /t, /m, /s. The post-/i
        //    layers are the *main* stereo layers.
        // 2. Isotope-specific stereo — isotopic substitution changes CIP priorities,
        //    producing a second /b, /t, /m, /s set after /i.
        //
        // If any sublayer type appears both before AND after /i, the
        // entire post-/i group is isotope-specific.  Otherwise the
        // post-/i layers fill in the (absent) main stereo.
        // Isotope-specific stereo may use shorthands (e.g. `m` for
        // mirror-image) that the parser does not yet fully support.
        // If parsing fails, restore the remainder so the content is
        // validated but not consumed.
        let saved_remainder = layer_remainder;
        let (pi_db, pi_tet, pi_allene, pi_info) = if let Ok(sublayers) =
            parse_stereo_sublayers(&mut layer_remainder, main_layer.chemical_formula())
        {
            sublayers
        } else {
            layer_remainder = saved_remainder;
            (None, None, None, None)
        };

        let has_duplicate = pi_db.is_some() && db.is_some()
            || pi_tet.is_some() && tet.is_some()
            || pi_allene.is_some() && allene.is_some()
            || pi_info.is_some() && s_info.is_some();

        let isotope_stereochemistry = if has_duplicate {
            build_stereo_layer(pi_db, pi_tet, pi_allene, pi_info)
        } else {
            if db.is_none() {
                db = pi_db;
            }
            if tet.is_none() {
                tet = pi_tet;
            }
            if allene.is_none() {
                allene = pi_allene;
            }
            if s_info.is_none() {
                s_info = pi_info;
            }
            None
        };

        let stereochemistry = build_stereo_layer(db, tet, allene, s_info);

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
            isotope_stereochemistry,
            fixed_hydrogen: None,
            reconnected: None,
            _version: core::marker::PhantomData,
        })
    }
}
