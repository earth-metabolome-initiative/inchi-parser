use crate::traits::prefix::Prefix;

pub struct IsotopeLayer;

impl Prefix for IsotopeLayer {
    const PREFIX: char = 'i';
}
