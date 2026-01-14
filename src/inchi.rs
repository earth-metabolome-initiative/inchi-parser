pub mod charge;
pub mod fixed_hydrogen;
pub mod isotope;
pub mod main_layer;
pub mod reconnected;
pub mod stereochemistry;
pub use charge::ChargeLayer;
pub use fixed_hydrogen::FixedHydrogenLayer;
pub use isotope::IsotopeLayer;
pub use main_layer::MainLayer;
pub use reconnected::ReconnectedLayer;
pub use stereochemistry::StereochemistryLayer;

use crate::version::Version;

pub struct InChI<V: Version = crate::version::StandardVersion1_07_1> {
    main_layer: MainLayer,
    charge: ChargeLayer,
    stereochemistry: StereochemistryLayer,
    isotope: IsotopeLayer,
    fixed_hydrogen: FixedHydrogenLayer,
    reconnected: ReconnectedLayer,
    _version: std::marker::PhantomData<V>,
}
