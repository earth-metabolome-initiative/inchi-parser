//! Module for the isotope layer of an InChI.

use alloc::vec::Vec;

use elements_rs::isotopes::HydrogenIsotope;

use crate::traits::prefix::Prefix;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// A non-hydrogen atom isotope specification.
pub struct IsotopeAtom {
    /// 0-based atom index.
    pub(crate) atom_index: u16,
    /// Mass shift relative to the rounded average atomic mass (can be
    /// negative).
    pub(crate) mass_shift: i16,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// A hydrogen isotope specification within the isotope layer.
pub struct IsotopeHydrogen {
    /// Which hydrogen isotope (D, T, or H1).
    pub(crate) isotope: HydrogenIsotope,
    /// How many of this isotope.
    pub(crate) count: u16,
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// A single component's isotope specifications.
pub struct IsotopeComponent {
    /// Non-hydrogen atom isotope specifications.
    pub(crate) atoms: Vec<IsotopeAtom>,
    /// Hydrogen isotope specifications.
    pub(crate) hydrogens: Vec<IsotopeHydrogen>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// The isotope layer of an InChI.
pub struct IsotopeLayer {
    pub(crate) components: Vec<IsotopeComponent>,
}

impl IsotopeAtom {
    /// Returns the 0-based atom index.
    #[must_use]
    pub fn atom_index(&self) -> u16 {
        self.atom_index
    }

    /// Returns the mass shift relative to the rounded average atomic mass.
    #[must_use]
    pub fn mass_shift(&self) -> i16 {
        self.mass_shift
    }
}

impl IsotopeHydrogen {
    /// Returns which hydrogen isotope (D, T, or H1).
    #[must_use]
    pub fn isotope(&self) -> HydrogenIsotope {
        self.isotope
    }

    /// Returns how many of this isotope.
    #[must_use]
    pub fn count(&self) -> u16 {
        self.count
    }
}

impl IsotopeComponent {
    /// Returns the non-hydrogen atom isotope specifications.
    #[must_use]
    pub fn atoms(&self) -> &[IsotopeAtom] {
        &self.atoms
    }

    /// Returns the hydrogen isotope specifications.
    #[must_use]
    pub fn hydrogens(&self) -> &[IsotopeHydrogen] {
        &self.hydrogens
    }
}

impl IsotopeLayer {
    /// Returns the per-component isotope specifications.
    #[must_use]
    pub fn components(&self) -> &[IsotopeComponent] {
        &self.components
    }
}

impl Prefix for IsotopeLayer {
    const PREFIX: char = 'i';
}
