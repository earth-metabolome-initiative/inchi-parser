/// Trait for parsing InChI layers.
pub trait ParseLayer: Sized {
    type Error;
    type Context<'a>;
    fn parse(input: &str, context: Self::Context<'_>) -> Result<Self, Self::Error>;
}
