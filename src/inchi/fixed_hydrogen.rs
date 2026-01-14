use crate::traits::prefix::Prefix;

pub struct FixedHydrogenLayer;

impl Prefix for FixedHydrogenLayer {
    const PREFIX: char = 'f';
}
