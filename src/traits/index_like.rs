//! Submodule defining traits for numbers that can be used as indices.

use core::ops::{AddAssign, SubAssign};

use geometric_traits::traits::IntoUsize;
use molecular_formulas::NumberLike;
use num_traits::{Bounded, SaturatingAdd, SaturatingSub, ToPrimitive, Unsigned};

/// A trait for types that can be used as indices.
pub trait IndexLike:
    NumberLike
    + Ord
    + IntoUsize
    + ToPrimitive
    + SaturatingAdd
    + SaturatingSub
    + Unsigned
    + Bounded
    + AddAssign
    + SubAssign
{
}

impl<T> IndexLike for T where
    T: NumberLike
        + Ord
        + IntoUsize
        + ToPrimitive
        + SaturatingAdd
        + SaturatingSub
        + Unsigned
        + Bounded
        + AddAssign
        + SubAssign
{
}
