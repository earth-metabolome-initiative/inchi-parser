//! Component-level parser for the hydrogen layer.

use alloc::vec::Vec;


use crate::{
    errors::HydrogenLayerTokenError,
    inchi::main_layer::{HydrogenComponent, MobileHydrogenGroup},
    traits::IndexLike,
};

use super::sub_tokens::{HydogenLayerSubTokens, HydrogenLayerSubTokenIter};

/// Parses one component string (after stripping any `n*` prefix) into a
/// [`HydrogenComponent`].
///
/// Atom indices are converted from 1-based to 0-based here.
/// `num_atoms` comes from the subformula's non-hydrogen atom count and is used
/// to pre-size the `fixed_h` array.
pub(super) fn parse_component<Idx>(
    input: &str,
    num_atoms: usize,
) -> Result<HydrogenComponent<Idx>, HydrogenLayerTokenError<Idx>>
where
    Idx: IndexLike,
    u8: TryFrom<Idx, Error = core::num::TryFromIntError>,
{
    let mut fixed_h: Vec<u8> = alloc::vec![0u8; num_atoms];
    let mut mobile_groups: Vec<MobileHydrogenGroup<Idx>> = Vec::new();

    // Pending atom indices (1-based) accumulated before the H token
    let mut atom_buf: Vec<Idx> = Vec::new();

    // Mobile-group state
    let mut in_mobile = false;
    let mut mobile_count = 0u8;
    let mut mobile_charged = false;
    let mut mobile_atoms: Vec<Idx> = Vec::new();

    let iter = HydrogenLayerSubTokenIter::<Idx>::from(input);

    for result in iter {
        let token = result?;
        if in_mobile {
            match token {
                HydogenLayerSubTokens::Comma => {
                    // Comma separates atom indices inside a mobile group — skip it
                }
                HydogenLayerSubTokens::Index(n) => {
                    mobile_atoms.push(n);
                }
                HydogenLayerSubTokens::CloseParenthesis => {
                    // Convert 1-based to 0-based and emit the group
                    let atoms = mobile_atoms
                        .drain(..)
                        .map(|a| a.saturating_sub(&Idx::ONE))
                        .collect();
                    mobile_groups.push(MobileHydrogenGroup {
                        count: mobile_count,
                        charged: mobile_charged,
                        atoms,
                    });
                    in_mobile = false;
                }
                HydogenLayerSubTokens::Range(_) => {
                    // Ranges are not valid inside mobile groups
                    return Err(HydrogenLayerTokenError::InvalidCharacter('-'));
                }
                HydogenLayerSubTokens::H(_) => {
                    return Err(HydrogenLayerTokenError::InvalidCharacter('H'));
                }
                HydogenLayerSubTokens::SharedHydrogens { .. } => {
                    return Err(HydrogenLayerTokenError::InvalidCharacter('('));
                }
                HydogenLayerSubTokens::Asterisk(_) => {
                    return Err(HydrogenLayerTokenError::InvalidCharacter('*'));
                }
            }
        } else {
            match token {
                HydogenLayerSubTokens::Comma => {
                    // Comma between fixed entries or between fixed and mobile blocks — skip it
                }
                HydogenLayerSubTokens::Index(n) => {
                    atom_buf.push(n);
                }
                HydogenLayerSubTokens::Range((s, e)) => {
                    // Expand the range into atom_buf (1-based indices)
                    let mut i = s;
                    loop {
                        atom_buf.push(i);
                        if i == e {
                            break;
                        }
                        i += Idx::ONE;
                    }
                }
                HydogenLayerSubTokens::H(count) => {
                    // Assign fixed H count to all pending atoms
                    for idx in atom_buf.drain(..) {
                        let zero_based = idx.saturating_sub(&Idx::ONE).into_usize();
                        if zero_based >= fixed_h.len() {
                            fixed_h.resize(zero_based + 1, 0);
                        }
                        fixed_h[zero_based] = count;
                    }
                }
                HydogenLayerSubTokens::SharedHydrogens { count, charged } => {
                    // Start a new mobile group
                    in_mobile = true;
                    mobile_count = count;
                    mobile_charged = charged;
                    mobile_atoms.clear();
                }
                HydogenLayerSubTokens::CloseParenthesis => {
                    return Err(HydrogenLayerTokenError::InvalidCharacter(')'));
                }
                HydogenLayerSubTokens::Asterisk(_) => {
                    // Asterisk should have been consumed by the caller before invoking
                    // parse_component
                    return Err(HydrogenLayerTokenError::InvalidCharacter('*'));
                }
            }
        }
    }

    // Unclosed mobile group
    if in_mobile {
        return Err(HydrogenLayerTokenError::UnexpectedEndOfInput(
            HydogenLayerSubTokens::CloseParenthesis,
        ));
    }

    // Pending atoms without a following H token
    if let Some(&last) = atom_buf.last() {
        return Err(HydrogenLayerTokenError::UnexpectedEndOfInput(
            HydogenLayerSubTokens::Index(last),
        ));
    }

    Ok(HydrogenComponent { fixed_h, mobile_groups })
}
