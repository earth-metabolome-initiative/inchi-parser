use core::str::FromStr;

use crate::{errors::Error, inchi::proton_layer::ProtonSublayer};

impl FromStr for ProtonSublayer {
    type Err = crate::errors::Error<usize>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Err(Error::UnimplementedFeature("Proton Sublayer not implemented yet !"))
    }
}
