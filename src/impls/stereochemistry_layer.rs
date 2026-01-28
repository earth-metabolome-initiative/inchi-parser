use core::str::FromStr;

use crate::{
    errors::Error,
    inchi::stereochemistry_layer::{
        AlleneSublayer, DoubleBondSublayer, StereoChemistryInformationSublayer,
        StereochemistryLayer, TetrahedralSublayer,
    },
};

impl FromStr for StereochemistryLayer {
    type Err = crate::errors::Error<usize>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Err(Error::UnimplementedFeature("Stereochemistry layer not implemented yet !"))
    }
}

impl FromStr for AlleneSublayer {
    type Err = crate::errors::Error<usize>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Err(Error::UnimplementedFeature("AlleneSublayer not implemented yet !"))
    }
}

impl FromStr for DoubleBondSublayer {
    type Err = crate::errors::Error<usize>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Err(Error::UnimplementedFeature("DoubleBondSublayer not implemented yet !"))
    }
}

impl FromStr for StereoChemistryInformationSublayer {
    type Err = crate::errors::Error<usize>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Err(Error::UnimplementedFeature("StereoChemistryInformationSublayer not implemented yet !"))
    }
}

impl FromStr for TetrahedralSublayer {
    type Err = crate::errors::Error<usize>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Err(Error::UnimplementedFeature("TetrahedralSublayer not implemented yet !"))
    }
}
