//! Module for the stereochemistry layer of an InChI.

use alloc::vec::Vec;

use crate::traits::prefix::Prefix;

/// Stereo parity designation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StereoParity {
    /// `+` parity (trans/E for double bonds, positive for tetrahedral).
    Plus,
    /// `-` parity (cis/Z for double bonds, negative for tetrahedral).
    Minus,
    /// `?` parity (unknown).
    Unknown,
}

/// A double bond stereo specification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DoubleBondStereo {
    pub(crate) atom1: u16,
    pub(crate) atom2: u16,
    pub(crate) parity: StereoParity,
}

impl DoubleBondStereo {
    /// Returns the first atom index (0-based).
    #[must_use]
    pub fn atom1(&self) -> u16 {
        self.atom1
    }

    /// Returns the second atom index (0-based).
    #[must_use]
    pub fn atom2(&self) -> u16 {
        self.atom2
    }

    /// Returns the parity.
    #[must_use]
    pub fn parity(&self) -> StereoParity {
        self.parity
    }
}

/// Double bond stereo sublayer: per-component bond specs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DoubleBondSublayer {
    pub(crate) components: Vec<Vec<DoubleBondStereo>>,
}

impl DoubleBondSublayer {
    /// Returns the per-component double bond stereo specifications.
    #[must_use]
    pub fn components(&self) -> &[Vec<DoubleBondStereo>] {
        &self.components
    }
}

impl Prefix for DoubleBondSublayer {
    const PREFIX: char = 'b';
}

/// A tetrahedral stereo specification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TetrahedralStereo {
    pub(crate) atom: u16,
    pub(crate) parity: StereoParity,
}

impl TetrahedralStereo {
    /// Returns the atom index (0-based).
    #[must_use]
    pub fn atom(&self) -> u16 {
        self.atom
    }

    /// Returns the parity.
    #[must_use]
    pub fn parity(&self) -> StereoParity {
        self.parity
    }
}

/// A single component's tetrahedral stereo content.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TetrahedralComponent {
    /// Explicit parity assignments for stereo centers.
    Explicit(Vec<TetrahedralStereo>),
    /// Abbreviation `m`: same stereo as the main (non-isotopic) layer.
    /// Only valid in isotope-specific layers.
    SameAsMainLayer,
}

/// Tetrahedral stereo sublayer: per-component center specs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TetrahedralSublayer {
    pub(crate) components: Vec<TetrahedralComponent>,
}

impl TetrahedralSublayer {
    /// Returns the per-component tetrahedral stereo specifications.
    #[must_use]
    pub fn components(&self) -> &[TetrahedralComponent] {
        &self.components
    }
}

impl Prefix for TetrahedralSublayer {
    const PREFIX: char = 't';
}

/// Allene/mirror sublayer: per-fragment parity groups.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlleneSublayer {
    pub(crate) values: Vec<Vec<u8>>,
}

impl AlleneSublayer {
    /// Returns the per-fragment parity groups.
    #[must_use]
    pub fn values(&self) -> &[Vec<u8>] {
        &self.values
    }
}

impl Prefix for AlleneSublayer {
    const PREFIX: char = 'm';
}

/// Stereo type sublayer: global value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoChemistryInformationSublayer {
    pub(crate) value: u8,
}

impl StereoChemistryInformationSublayer {
    /// Returns the stereo type value (1 = absolute, 2 = relative, 3 = racemic).
    #[must_use]
    pub fn value(&self) -> u8 {
        self.value
    }
}

impl Prefix for StereoChemistryInformationSublayer {
    const PREFIX: char = 's';
}

/// Composite stereochemistry layer holding all sublayers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereochemistryLayer {
    pub(crate) double_bond: Option<DoubleBondSublayer>,
    pub(crate) tetrahedral: Option<TetrahedralSublayer>,
    pub(crate) allene: Option<AlleneSublayer>,
    pub(crate) stereo_info: Option<StereoChemistryInformationSublayer>,
}

impl StereochemistryLayer {
    /// Returns the double bond stereo sublayer, if present.
    #[must_use]
    pub fn double_bond(&self) -> Option<&DoubleBondSublayer> {
        self.double_bond.as_ref()
    }

    /// Returns the tetrahedral stereo sublayer, if present.
    #[must_use]
    pub fn tetrahedral(&self) -> Option<&TetrahedralSublayer> {
        self.tetrahedral.as_ref()
    }

    /// Returns the allene/mirror sublayer, if present.
    #[must_use]
    pub fn allene(&self) -> Option<&AlleneSublayer> {
        self.allene.as_ref()
    }

    /// Returns the stereo type sublayer, if present.
    #[must_use]
    pub fn stereo_info(&self) -> Option<&StereoChemistryInformationSublayer> {
        self.stereo_info.as_ref()
    }
}
