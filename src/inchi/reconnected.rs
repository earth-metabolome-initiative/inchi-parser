use crate::traits::prefix::Prefix;

pub struct ReconnectedLayer;

impl Prefix for ReconnectedLayer {
    const PREFIX: char = 'r';
}
