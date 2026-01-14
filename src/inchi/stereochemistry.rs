use crate::traits::prefix::Prefix;

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
