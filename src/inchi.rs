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

use crate::{
    inchi::{charge_layer::ChargeSubLayer, proton_layer::ProtonSublayer},
    version::Version,
};

impl<V: Version> InChI<V> {
    /// Returns the per-component charges from the `/q` layer, if present.
    #[must_use]
    pub fn charges(&self) -> Option<&[i16]> {
        self.charge.as_ref().map(|c| c.charges.as_slice())
    }

    /// Returns the global proton balance from the `/p` layer, if present.
    #[must_use]
    pub fn proton_count(&self) -> Option<i16> {
        self.proton.as_ref().map(|p| p.proton_count)
    }

    /// Returns the isotope layer, if present.
    #[must_use]
    pub fn isotope(&self) -> Option<&IsotopeLayer> {
        self.isotope.as_ref()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// The InChI structure
pub struct InChI<V: Version = crate::version::StandardVersion1_07_4> {
    pub(crate) main_layer: MainLayer,
    pub(crate) charge: Option<ChargeSubLayer>,
    pub(crate) proton: Option<ProtonSublayer>,
    pub(crate) stereochemistry: Option<StereochemistryLayer>,
    pub(crate) isotope: Option<IsotopeLayer>,
    pub(crate) fixed_hydrogen: Option<FixedHydrogenLayer>,
    pub(crate) reconnected: Option<ReconnectedLayer>,
    pub(crate) _version: core::marker::PhantomData<V>,
}
