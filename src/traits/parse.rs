use crate::traits::prefix::Prefix;

/// Trait for parsing InChI layers.

pub trait FromStrWithContext: Sized {
    type Context<'a>;
    /// Given a string that matches a specific layer, the function tries to
    /// parse that layer and returns the struct for that layer or an error.
    fn from_str_with_context(
        input: &str,
        context: Self::Context<'_>,
    ) -> Result<Self, crate::errors::Error<usize>>;
}

pub trait ConsumeStr: Sized {
    /// Takes a string as input and checks if it matches the layer.
    /// If it does, then the begining of the string that matches a layer is
    /// consumed. If it does then the string return itself.
    fn consume_str(input: &str) -> Result<(Self, &str), crate::errors::Error<usize>>;
}

pub trait PrefixFromStrWithContext: Prefix + FromStrWithContext {
    fn try_build_layer(
        input: &mut &str,
        context: Self::Context<'_>,
    ) -> Result<Option<Self>, crate::errors::Error<usize>> {
        let layer = if input.starts_with(Self::PREFIX) {
            let (layer_string, layer_remainder_string) =
                input.split_once('/').unwrap_or((input, ""));
            *input = layer_remainder_string;
            Some(Self::from_str_with_context(layer_string, context)?)
        } else {
            None
        };
        Ok(layer)
    }
}
