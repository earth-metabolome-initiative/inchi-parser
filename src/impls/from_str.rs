use crate::inchi::InChI;
use crate::version::Version;
use std::str::FromStr;

impl<V: Version> FromStr for InChI<V> {
    type Err = crate::errors::Errors;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // First we remove the "InChI=" prefix
        let Some(s) = s.strip_prefix(crate::constants::INCHI_PREFIX) else {
            return Err(crate::errors::Errors::MissingInchiPrefix);
        };

        // Next the version prefix is removed
        let Some(s) = s.strip_prefix(V::VERSION_PREFIX) else {
            return Err(crate::errors::Errors::MissingVersionPrefix);
        };

        // Then we remove the first '/' to get the first layer
        let Some(s) = s.strip_prefix('/') else {
            return Err(crate::errors::Errors::MissingForwardSlash);
        };
    }
}
