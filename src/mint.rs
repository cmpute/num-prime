//! Wrapper of integer to makes it efficient in modular arithmetics but still have the same
//! API of normal integers.

use core::ops::*;
use either::*;
use num_integer::Integer;
use num_modular::{ModularInteger, Montgomery, MontgomeryInt};

// TODO (v0.3.x): Implement basic support for these

/// Integer with fast modular arithmetics support, based on MontgomeryInt
pub struct Mint<T: Integer + Montgomery>(Either<T, MontgomeryInt<T>>);

impl<T: Integer + Montgomery> From<T> for Mint<T> {
    fn from(v: T) -> Self {
        Self(Left(v))
    }
}

impl<T: Integer + Montgomery + Clone> Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    pub fn value(&self) -> T {
        match &self.0 {
            Left(v) => v.clone(),
            Right(m) => m.residue(),
        }
    }
}

impl<T: Integer + Montgomery + Clone> Add for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(match (self.0, rhs.0) {
            (Left(v1), Left(v2)) => Left(v1.add(v2)),
            (Left(v1), Right(v2)) => Right(v2.convert(v1).add(v2)),
            (Right(v1), Left(v2)) => {
                let v2 = v1.convert(v2);
                Right(v1.add(v2))
            }
            (Right(v1), Right(v2)) => Right(v1.add(v2)),
        })
    }
}
