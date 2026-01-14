/// Trait for parsing InChI layers.
pub trait ParseLayer: Sized {
    type Error;
    type Context;
    fn parse(input: &str, context: &mut Self::Context) -> Result<Self, Self::Error>;
}
