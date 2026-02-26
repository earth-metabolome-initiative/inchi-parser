//! Module for the hydrogen sub-layer of an InChI main layer.

use alloc::vec::Vec;

use crate::traits::prefix::Prefix;

/// A group of mobile (tautomeric) hydrogens delocalized over `atoms`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct MobileHydrogenGroup<Idx = u16> {
    /// Number of mobile hydrogens in this group.
    pub(crate) count: u8,
    /// Whether these are charged mobile hydrogens (H⁻, from `(H-,...)`).
    pub(crate) charged: bool,
    /// Zero-based atom indices that are possible protonation sites.
    pub(crate) atoms: Vec<Idx>,
}

/// Hydrogen data for one molecular component.
///
/// `fixed_h[i]` is indexed by **0-based atom index** within this component.
/// The value is the number of fixed hydrogens on atom `i`.
/// `0` means "not mentioned in `/h`" — atoms with 0 H are absent from the layer.
///
/// The value type is `u8` because it is a **hydrogen count** (always 0–8),
/// not an atom index.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct HydrogenComponent<Idx = u16> {
    /// Per-atom fixed hydrogen counts, indexed by 0-based atom index.
    /// `0` means no fixed hydrogens.
    pub(crate) fixed_h: Vec<u8>,
    /// Mobile hydrogen groups for this component.
    pub(crate) mobile_groups: Vec<MobileHydrogenGroup<Idx>>,
}

/// The parsed hydrogen sublayer: one `HydrogenComponent` per `;`-delimited fragment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct HydrogensSubLayer {
    /// Per-fragment hydrogen data, in the same order as formula/connection layer components.
    pub(crate) components: Vec<HydrogenComponent>,
}

impl Prefix for HydrogensSubLayer {
    const PREFIX: char = 'h';
}
