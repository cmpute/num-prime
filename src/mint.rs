//! Wrapper of integer to makes it efficient in modular arithmetics but still have the same
//! API of normal integers.

use core::ops::*;
use either::*;
use num_traits::{Num, Zero, One};
use num_integer::Integer;
use num_modular::{ModularInteger, Montgomery, MontgomeryInt};

// TODO (v0.3.x): Implement PrimalityBase for Mint

/// Integer with fast modular arithmetics support, based on MontgomeryInt
pub struct Mint<T: Integer + Montgomery>(Either<T, MontgomeryInt<T>>);

// it seems that auto derivation struggles to provide an implementation for Copy, Clone and Debug with proper trait bounds
impl<T: Integer + Montgomery + Clone> Clone for Mint<T> where T::Inv : Clone {
    #[inline(always)]
    fn clone(&self) -> Self {
        Self(self.0.clone())
    }
}
impl<T: Integer + Montgomery + Copy> Copy for Mint<T> where T::Inv : Copy { }
impl<T: Integer + Montgomery + core::fmt::Debug> core::fmt::Debug for Mint<T> where T::Inv : core::fmt::Debug {
    #[inline(always)]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl<T: Integer + Montgomery> From<T> for Mint<T> {
    #[inline(always)]
    fn from(v: T) -> Self {
        Self(Left(v))
    }
}

#[inline(always)]
fn left_only<T: Integer + Montgomery>(lhs: Mint<T>, rhs: Mint<T>) -> (T, T) {
    match(lhs.0, rhs.0) {
        (Left(v1), Left(v2)) => (v1, v2),
        (_, _) => panic!("not supported"),
    }
}

#[inline(always)]
fn left_ref_only<'a, T: Integer + Montgomery>(lhs: &'a Mint<T>, rhs: &'a Mint<T>) -> (&'a T, &'a T) {
    match(&lhs.0, &rhs.0) {
        (Left(v1), Left(v2)) => (v1, v2),
        (_, _) => panic!("not supported"),
    }
}

macro_rules! forward_binops_left_ref_only {
    ($method:ident) => {
        #[inline(always)]
        fn $method(&self, other: &Self) -> Self {
            let (v1, v2) = left_ref_only(self, other);
            Self(Left(v1.$method(v2)))
        }
    };
    ($method:ident => $return:ty) => {
        #[inline(always)]
        fn $method(&self, other: &Self) -> $return {
            let (v1, v2) = left_ref_only(self, other);
            v1.$method(v2)
        }
    };
}

macro_rules! forward_uops_ref {
    ($method:ident => $return:ty) => {
        #[inline(always)]
        fn $method(&self) -> $return {
            match &self.0 {
                Left(v) => v.$method(),
                Right(m) => m.residue().$method()
            }
        }
    };
}

impl<T: Integer + Montgomery + Clone> PartialEq for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone, {
    fn eq(&self, other: &Self) -> bool {
        match (&self.0, &other.0) {
            (Left(v1), Left(v2)) => v1 == v2,
            (Left(v1), Right(v2)) => v1 == &v2.residue(),
            (Right(v1), Left(v2)) => &v1.residue() == v2,
            (Right(v1), Right(v2)) => {
                debug_assert!(v1.modulus() == v2.modulus());
                v1.residue() == v2.residue()
            },
        }
    }
}
impl<T: Integer + Montgomery + Clone> Eq for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone, { }

impl<T: Integer + Montgomery + Clone> PartialOrd for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone, {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match (&self.0, &other.0) {
            (Left(v1), Left(v2)) => v1.partial_cmp(v2),
            (Left(v1), Right(v2)) => v1.partial_cmp(&v2.residue()),
            (Right(v1), Left(v2)) => v1.residue().partial_cmp(v2),
            (Right(v1), Right(v2)) => {
                debug_assert!(v1.modulus() == v2.modulus());
                v1.residue().partial_cmp(&v2.residue())
            },
        }
    }
}
impl<T: Integer + Montgomery + Clone> Ord for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{ 
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (&self.0, &other.0) {
            (Left(v1), Left(v2)) => v1.cmp(v2),
            (Left(v1), Right(v2)) => v1.cmp(&v2.residue()),
            (Right(v1), Left(v2)) => v1.residue().cmp(v2),
            (Right(v1), Right(v2)) => v1.residue().cmp(&v2.residue()),
        }
    }
}

impl<T: Integer + Montgomery + Clone> Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    #[inline(always)]
    pub fn value(&self) -> T {
        match &self.0 {
            Left(v) => v.clone(),
            Right(m) => m.residue(),
        }
    }
}

// forward binary operators by converting result to MontgomeryInt whenever possible
macro_rules! forward_binops_right {
    (impl $imp:ident, $method:ident) => {
        impl<T: Integer + Montgomery + Clone> $imp for Mint<T>
        where
            T::Double: From<T>,
            T::Inv: Clone,
        {
            type Output = Self;

            #[inline]
            fn $method(self, rhs: Self) -> Self::Output {
                Self(match (self.0, rhs.0) {
                    (Left(v1), Left(v2)) => Left(v1.$method(v2)),
                    (Left(v1), Right(v2)) => Right(v2.convert(v1).$method(v2)),
                    (Right(v1), Left(v2)) => {
                        let v2 = v1.convert(v2);
                        Right(v1.$method(v2))
                    }
                    (Right(v1), Right(v2)) => Right(v1.$method(v2)),
                })
            }
        }
    }
}

forward_binops_right!(impl Add, add);
forward_binops_right!(impl Sub, sub);
forward_binops_right!(impl Mul, mul);

impl<T: Integer + Montgomery + Clone> Div for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        let (v1, v2) = left_only(self, rhs);
        Self(Left(v1.div(v2)))
    }
}

impl<T: Integer + Montgomery + Clone> Rem for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    type Output = Self;

    #[inline]
    fn rem(self, rhs: Self) -> Self::Output {
        let (v1, v2) = left_only(self, rhs);
        Self(Right(MontgomeryInt::new(v1, v2)))
    }
}

impl<T: Integer + Montgomery + Clone> Zero for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    #[inline(always)]
    fn zero() -> Self {
        Self(Left(T::zero()))
    }
    #[inline(always)]
    fn is_zero(&self) -> bool {
        match &self.0 {
            Left(v) => v.is_zero(),
            Right(m) => m.residue().is_zero() // TODO(v0.3.1): Add is_zero for MontgomeryInt
        }
    }
}

impl<T: Integer + Montgomery + Clone> One for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    #[inline(always)]
    fn one() -> Self {
        Self(Left(T::one()))
    }
    forward_uops_ref!(is_one => bool);
}

impl<T: Integer + Montgomery + Clone> Num for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    type FromStrRadixErr = <T as Num>::FromStrRadixErr;

    #[inline(always)]
    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        T::from_str_radix(str, radix).map(|v| Self(Left(v)))
    }
}

impl<T: Integer + Montgomery + Clone> Integer for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    forward_binops_left_ref_only!(div_floor);
    forward_binops_left_ref_only!(mod_floor);
    forward_binops_left_ref_only!(gcd);
    forward_binops_left_ref_only!(lcm);
    forward_binops_left_ref_only!(divides => bool);
    forward_binops_left_ref_only!(is_multiple_of => bool);
    forward_uops_ref!(is_even => bool);
    forward_uops_ref!(is_odd => bool);

    #[inline(always)]
    fn div_rem(&self, other: &Self) -> (Self, Self) {
        let (v1, v2) = left_ref_only(self, other);
        let (q, r) = v1.div_rem(v2);
        (Self(Left(q)), Self(Left(r)))
    }
}
