//! Module for the stereochemistry layer of an InChI.
use crate::traits::prefix::Prefix;

#[derive(Debug, Clone, PartialEq, Eq)]
/// The stereochemistry layer of an InChI.
pub(crate) struct StereochemistryLayer;

#[expect(dead_code, reason = "stereo sublayer parsing is not implemented yet")]
pub(crate) struct DoubleBondSublayer;
#[expect(dead_code, reason = "stereo sublayer parsing is not implemented yet")]
pub(crate) struct TetrahedralSublayer;
#[expect(dead_code, reason = "stereo sublayer parsing is not implemented yet")]
pub(crate) struct AlleneSublayer;
#[expect(dead_code, reason = "stereo sublayer parsing is not implemented yet")]
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
