//! Component-level parser for the hydrogen layer.

use alloc::vec::Vec;

use crate::{
    errors::HydrogenLayerTokenError,
    inchi::main_layer::{HydrogenComponent, MobileHydrogenGroup},
    traits::IndexLike,
};

use super::sub_tokens::{HydrogenLayerSubTokens, HydrogenLayerSubTokenIter};

/// Validates a 1-based atom index and converts it to 0-based, checking both
/// for zero (InChI indices start at 1) and for exceeding `num_atoms`.
fn validate_atom_index<Idx: IndexLike>(
    idx: Idx,
    num_atoms: usize,
) -> Result<usize, HydrogenLayerTokenError<Idx>> {
    if idx < Idx::ONE {
        return Err(HydrogenLayerTokenError::ZeroAtomIndex);
    }
    let zero_based = (idx - Idx::ONE).into_usize();
    if zero_based >= num_atoms {
        return Err(HydrogenLayerTokenError::AtomIndexOutOfBounds {
            index: idx,
            num_atoms,
        });
    }
    Ok(zero_based)
}

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
                HydrogenLayerSubTokens::Comma => {
                    // Comma separates atom indices inside a mobile group — skip it
                }
                HydrogenLayerSubTokens::Index(n) => {
                    validate_atom_index(n, num_atoms)?;
                    mobile_atoms.push(n);
                }
                HydrogenLayerSubTokens::CloseParenthesis => {
                    // Convert 1-based to 0-based and emit the group
                    // (indices already validated above)
                    let atoms = mobile_atoms
                        .drain(..)
                        .map(|a| a - Idx::ONE)
                        .collect();
                    mobile_groups.push(MobileHydrogenGroup {
                        count: mobile_count,
                        charged: mobile_charged,
                        atoms,
                    });
                    in_mobile = false;
                }
                HydrogenLayerSubTokens::Range(_) => {
                    // Ranges are not valid inside mobile groups
                    return Err(HydrogenLayerTokenError::InvalidCharacter('-'));
                }
                HydrogenLayerSubTokens::H(_) => {
                    return Err(HydrogenLayerTokenError::InvalidCharacter('H'));
                }
                HydrogenLayerSubTokens::SharedHydrogens { .. } => {
                    return Err(HydrogenLayerTokenError::InvalidCharacter('('));
                }
                HydrogenLayerSubTokens::Asterisk(_) => {
                    return Err(HydrogenLayerTokenError::InvalidCharacter('*'));
                }
            }
        } else {
            match token {
                HydrogenLayerSubTokens::Comma => {
                    // Comma between fixed entries or between fixed and mobile blocks — skip it
                }
                HydrogenLayerSubTokens::Index(n) => {
                    atom_buf.push(n);
                }
                HydrogenLayerSubTokens::Range((s, e)) => {
                    if s > e {
                        return Err(HydrogenLayerTokenError::InvalidRange(s, e));
                    }
                    // Validate both endpoints
                    validate_atom_index(s, num_atoms)?;
                    validate_atom_index(e, num_atoms)?;
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
                HydrogenLayerSubTokens::H(count) => {
                    // Assign fixed H count to all pending atoms
                    for idx in atom_buf.drain(..) {
                        let zero_based = validate_atom_index(idx, num_atoms)?;
                        fixed_h[zero_based] = count;
                    }
                }
                HydrogenLayerSubTokens::SharedHydrogens { count, charged } => {
                    // Start a new mobile group
                    in_mobile = true;
                    mobile_count = count;
                    mobile_charged = charged;
                    mobile_atoms.clear();
                }
                HydrogenLayerSubTokens::CloseParenthesis => {
                    return Err(HydrogenLayerTokenError::InvalidCharacter(')'));
                }
                HydrogenLayerSubTokens::Asterisk(_) => {
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
            HydrogenLayerSubTokens::CloseParenthesis,
        ));
    }

    // Pending atoms without a following H token
    if let Some(&last) = atom_buf.last() {
        return Err(HydrogenLayerTokenError::UnexpectedEndOfInput(
            HydrogenLayerSubTokens::Index(last),
        ));
    }

    Ok(HydrogenComponent { fixed_h, mobile_groups })
}
