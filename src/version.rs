/// Trait for the different InChI versions.
pub trait Version: Default {
    const VERSION: &'static str;
    const VERSION_PREFIX: &'static str;
}

/// InChI version 1.07.4
#[derive(Default)]
pub struct StandardVersion1_07_4;

#[derive(Default)]
pub struct Version1_07_4;

impl Version for StandardVersion1_07_4 {
    const VERSION: &'static str = "1.07.4";
    const VERSION_PREFIX: &'static str = "1S";
}

impl Version for Version1_07_4 {
    const VERSION: &'static str = "1.07.4";
    const VERSION_PREFIX: &'static str = "1";
}
