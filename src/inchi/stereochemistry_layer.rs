//! Module for the stereochemistry layer of an InChI.
use crate::traits::prefix::Prefix;

#[derive(Debug, Clone, PartialEq, Eq)]
/// The stereochemistry layer of an InChI.
pub(crate) struct StereochemistryLayer;

pub(crate) struct DoubleBondSublayer;
pub(crate) struct TetrahedralSublayer;
pub(crate) struct AlleneSublayer;
pub(crate) struct StereoChemistryInformationSublayer;

impl Prefix for DoubleBondSublayer {
    const PREFIX: char = 'b';
}

impl Prefix for TetrahedralSublayer {
    const PREFIX: char = 't';
}

impl Prefix for AlleneSublayer {
    const PREFIX: char = 'm';
}

impl Prefix for StereoChemistryInformationSublayer {
    const PREFIX: char = 's';
}
