//! Submodule defining traits for numbers that can be used as indices.

use core::{
    fmt::{Debug, Display},
    ops::Sub,
    str::FromStr,
};

use num_traits::{ConstOne, ConstZero};

/// A trait for types that can be used as indices.
pub trait IndexLike:
    Copy
    + Clone
    + Debug
    + Display
    + ConstZero
    + ConstOne
    + Sub<Output = Self>
    + PartialEq
    + Eq
    + FromStr
    + PartialOrd
    + Ord
{
}

impl<T> IndexLike for T where
    T: Copy
        + Clone
        + Debug
        + Display
        + ConstZero
        + ConstOne
        + Sub<Output = Self>
        + PartialEq
        + Eq
        + FromStr
        + PartialOrd
        + Ord
{
}
