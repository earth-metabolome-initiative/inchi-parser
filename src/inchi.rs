//! Module for the InChI structure and its layers.

pub mod charge_layer;
pub mod fixed_hydrogen;
pub mod isotope_layer;
pub mod main_layer;
pub mod proton_layer;
pub mod reconnected_layer;
pub mod stereochemistry_layer;
pub use fixed_hydrogen::FixedHydrogenLayer;
pub use isotope_layer::IsotopeLayer;
pub use main_layer::MainLayer;
pub use reconnected_layer::ReconnectedLayer;
pub(crate) use stereochemistry_layer::StereochemistryLayer;

use crate::{inchi::charge_layer::ChargeSubLayer, version::Version};

#[derive(Debug, Clone, PartialEq, Eq)]
/// The InChI structure
pub struct InChI<V: Version = crate::version::StandardVersion1_07_4> {
    main_layer: MainLayer,
    charge: ChargeSubLayer,
    stereochemistry: StereochemistryLayer,
    isotope: IsotopeLayer,
    fixed_hydrogen: FixedHydrogenLayer,
    reconnected: ReconnectedLayer,
    _version: core::marker::PhantomData<V>,
}
