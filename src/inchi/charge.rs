use crate::traits::prefix::Prefix;
pub struct ChargeLayer {
    charged_sublayer: Option<ChargeSubLayer>,
    proton_sublayer: Option<ProtonSublayer>,
}

struct ChargeSubLayer;
struct ProtonSublayer;

impl Prefix for ChargeSubLayer {
    const PREFIX: char = 'q';
}

impl Prefix for ProtonSublayer {
    const PREFIX: char = 'p';
}
