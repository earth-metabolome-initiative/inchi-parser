use core::str::Chars;

use crate::traits::{IndexLike, prefix::Prefix};

/// Trait for parsing InChI layers.
pub trait FromStrWithContext: Sized {
    type Idx: IndexLike;
    type Context<'a>;
    type Input<'a>: FromInChIStr<'a>;
    /// Given a string that matches a specific layer, the function tries to
    /// parse that layer and returns the struct for that layer or an error.
    fn from_str_with_context(
        input: Self::Input<'_>,
        context: Self::Context<'_>,
    ) -> Result<Self, crate::errors::Error<Self::Idx>>;
}

pub trait ConsumeStr: Sized {
    type Idx: IndexLike;
    /// Takes a string as input and checks if it matches the layer.
    /// If it does, then the begining of the string that matches a layer is
    /// consumed. If it does then the string return itself.
    fn consume_str(input: &str) -> Result<(Self, &str), crate::errors::Error<Self::Idx>>;
}

pub trait PrefixFromStrWithContext: Prefix + FromStrWithContext {
    fn try_build_layer(
        input: &mut &str,
        context: Self::Context<'_>,
    ) -> Result<Option<Self>, crate::errors::Error<Self::Idx>> {
        let layer = if input.starts_with(Self::PREFIX) {
            let (layer_string, layer_remainder_string) =
                input.split_once('/').unwrap_or((input, ""));
            *input = layer_remainder_string;
            Some(Self::from_str_with_context(FromInChIStr::from_inchi_str(layer_string), context)?)
        } else {
            None
        };
        Ok(layer)
    }
}

pub trait FromInChIStr<'a> {
    fn from_inchi_str(input: &'a str) -> Self;
}

impl<'a> FromInChIStr<'a> for &'a str {
    fn from_inchi_str(input: &'a str) -> Self {
        input
    }
}

impl<'a> FromInChIStr<'a> for core::iter::Peekable<Chars<'a>> {
    fn from_inchi_str(input: &'a str) -> Self {
        input.chars().peekable()
    }
}
