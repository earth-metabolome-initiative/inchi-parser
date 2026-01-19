//! Module for the stereochemistry layer of an InChI.
use crate::traits::prefix::Prefix;

#[derive(Debug, Clone, PartialEq, Eq)]
/// The stereochemistry layer of an InChI.
pub struct StereochemistryLayer;

struct DoubleBondSublayer;
struct TetrahedralSublayer;
struct AlleneSublayer;
struct StereoChemistryInformationSublayer;

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
