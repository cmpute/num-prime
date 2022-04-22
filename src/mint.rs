//! Wrapper of integer to makes it efficient in modular arithmetics but still have the same
//! API of normal integers.

use core::{ops::*, panic};
use either::*;
use num_traits::{Num, Zero, One, FromPrimitive, ToPrimitive, Pow};
use num_integer::{Integer, Roots};
use num_modular::{ModularInteger, Montgomery, MontgomeryInt};

use crate::{ExactRoots, BitTest};

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

        impl<T: Integer + Montgomery + Clone + for<'r> $imp<&'r T, Output = T>> $imp<&Self> for Mint<T>
        where
            T::Double: From<T>,
            T::Inv: Clone,
        {
            type Output = Mint<T>;
            #[inline]
            fn $method(self, rhs: &Self) -> Self::Output {
                Mint(match (self.0, &rhs.0) {
                    (Left(v1), Left(v2)) => Left(v1.$method(v2)),
                    (Left(v1), Right(v2)) => Right(v2.convert(v1).$method(v2)),
                    (Right(v1), Left(v2)) => {
                        let v2 = v1.convert(v2.clone());
                        Right(v1.$method(v2))
                    }
                    (Right(v1), Right(v2)) => Right(v1.$method(v2)),
                })
            }
        }

        impl<T: Integer + Montgomery + Clone> $imp<Mint<T>> for &Mint<T>
        where
            T::Double: From<T>,
            T::Inv: Clone,
        {
            type Output = Mint<T>;
            // FIXME: additional clone here due to https://github.com/rust-lang/rust/issues/39959
            // (same for ref & ref operation below, and those for Div and Rem)
            #[inline]
            fn $method(self, rhs: Mint<T>) -> Self::Output {
                Mint(match (&self.0, rhs.0) {
                    (Left(v1), Left(v2)) => Left(v1.clone().$method(v2)),
                    (Left(v1), Right(v2)) => Right(v2.convert(v1.clone()).$method(v2)),
                    (Right(v1), Left(v2)) => {
                        let v2 = v1.convert(v2);
                        Right(v1.clone().$method(v2))
                    }
                    (Right(v1), Right(v2)) => Right(v1.$method(v2)),
                })
            }
        }
        impl<T: Integer + Montgomery + Clone + for<'r> $imp<&'r T, Output = T>> $imp<Self> for &Mint<T>
        where
            T::Double: From<T>,
            T::Inv: Clone,
        {
            type Output = Mint<T>;
            #[inline]
            fn $method(self, rhs: Self) -> Self::Output {
                Mint(match (&self.0, &rhs.0) {
                    (Left(v1), Left(v2)) => Left(v1.clone().$method(v2)),
                    (Left(v1), Right(v2)) => Right(v2.convert(v1.clone()).$method(v2)),
                    (Right(v1), Left(v2)) => {
                        let v2 = v1.convert(v2.clone());
                        Right(v1.clone().$method(v2))
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
impl<T: Integer + Montgomery + Clone + for<'r> Div<&'r T, Output=T>> Div<&Self> for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    type Output = Self;

    #[inline]
    fn div(self, rhs: &Self) -> Self::Output {
        match(self.0, &rhs.0) {
            (Left(v1), Left(v2)) => Self(Left(v1.div(v2))),
            (_, _) => panic!("not supported"),
        }
    }
}
impl<T: Integer + Montgomery + Clone> Div<Mint<T>> for &Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    type Output = Mint<T>;

    #[inline]
    fn div(self, rhs: Mint<T>) -> Self::Output {
        match(&self.0, rhs.0) {
            (Left(v1), Left(v2)) => Mint(Left(v1.clone().div(v2))),
            (_, _) => panic!("not supported"),
        }
    }
}
impl<T: Integer + Montgomery + Clone + for<'r> Div<&'r T, Output=T>> Div<Self> for &Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    type Output = Mint<T>;

    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        match(&self.0, &rhs.0) {
            (Left(v1), Left(v2)) => Mint(Left(v1.clone().div(v2))),
            (_, _) => panic!("not supported"),
        }
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
        match(self.0, rhs.0) {
            (Left(v1), Left(v2)) => Self(Right(MontgomeryInt::new(v1, v2))),
            (Right(v1), Left(v2)) => {
                debug_assert!(v1.modulus() == &v2);
                Self(Right(v1))
            }
            (_, _) => panic!("not supported"),
        }
    }
}
impl<T: Integer + Montgomery + Clone + for<'r> Rem<&'r T, Output=T>> Rem<&Self> for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    type Output = Self;

    #[inline]
    fn rem(self, rhs: &Self) -> Self::Output {
        match(self.0, &rhs.0) {
            (Left(v1), Left(v2)) => Self(Right(MontgomeryInt::new(v1, v2.clone()))),
            (Right(v1), Left(v2)) => {
                debug_assert!(v1.modulus() == v2);
                Self(Right(v1))
            }
            (_, _) => panic!("not supported"),
        }
    }
}
impl<T: Integer + Montgomery + Clone> Rem<Mint<T>> for &Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    type Output = Mint<T>;

    #[inline]
    fn rem(self, rhs: Mint<T>) -> Self::Output {
        match(&self.0, rhs.0) {
            (Left(v1), Left(v2)) => Mint(Right(MontgomeryInt::new(v1.clone(), v2))),
            (Right(v1), Left(v2)) => {
                debug_assert!(v1.modulus() == &v2);
                Mint(Right(v1.clone()))
            }
            (_, _) => panic!("not supported"),
        }
    }
}
impl<T: Integer + Montgomery + Clone + for<'r> Rem<&'r T, Output=T>> Rem<Self> for &Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    type Output = Mint<T>;

    #[inline]
    fn rem(self, rhs: Self) -> Self::Output {
        match(&self.0, &rhs.0) {
            (Left(v1), Left(v2)) => Mint(Right(MontgomeryInt::new(v1.clone(), v2.clone()))),
            (Right(v1), Left(v2)) => {
                debug_assert!(v1.modulus() == v2);
                Mint(Right(v1.clone()))
            }
            (_, _) => panic!("not supported"),
        }
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
            Right(m) => m.is_zero()
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

impl<T: Integer + Montgomery + Clone + Roots> Roots for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone, {
    fn nth_root(&self, n: u32) -> Self {
        match &self.0 {
            Left(v) => Self(Left(v.nth_root(n))),
            Right(_) => panic!("not supported")
        }
    }
}

impl<T: Integer + Montgomery + Clone + FromPrimitive> FromPrimitive for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone, {
        fn from_f64(n: f64) -> Option<Self> {
            T::from_f64(n).map(|v| Self(Left(v)))
        }
        fn from_i64(n: i64) -> Option<Self> {
            T::from_i64(n).map(|v| Self(Left(v)))
        }
        fn from_u64(n: u64) -> Option<Self> {
            T::from_u64(n).map(|v| Self(Left(v)))
        }
}

impl<T: Integer + Montgomery + Clone + ToPrimitive> ToPrimitive for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone, {
        fn to_f64(&self) -> Option<f64> {
            match &self.0 {
                Left(v) => v.to_f64(),
                Right(m) => m.residue().to_f64()
            }
        }
        fn to_i64(&self) -> Option<i64> {
            match &self.0 {
                Left(v) => v.to_i64(),
                Right(m) => m.residue().to_i64()
            }
        }
        fn to_u64(&self) -> Option<u64> {
            match &self.0 {
                Left(v) => v.to_u64(),
                Right(m) => m.residue().to_u64()
            }
        }
}

impl<T: Integer + Montgomery + Clone + Pow<u32, Output = T>> Pow<u32> for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone, {
        type Output = Self;
        fn pow(self, rhs: u32) -> Self::Output {
            match self.0 {
                Left(v) => Self(Left(v.pow(rhs))),
                Right(_) => panic!("not supported")
            }
        }
}

impl<T: Integer + Montgomery + Clone + ExactRoots> ExactRoots for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone, {
        fn nth_root_exact(&self, n: u32) -> Option<Self> {
            match &self.0 {
                Left(v) => v.nth_root_exact(n).map(|v| Self(Left(v))),
                Right(_) => panic!("not supported")
            }
        }
}

impl<T: Integer + Montgomery + Clone + BitTest> BitTest for Mint<T> 
where
    T::Double: From<T>,
    T::Inv: Clone, {
   fn bit(&self, position: usize) -> bool {
        match &self.0 {
            Left(v) => v.bit(position),
            Right(_) => panic!("not supported")
        }
   } 
   fn bits(&self) -> usize {
        match &self.0 {
            Left(v) => v.bits(),
            Right(_) => panic!("not supported")
        }
   }
   fn trailing_zeros(&self) -> usize {
        match &self.0 {
            Left(v) => v.trailing_zeros(),
            Right(_) => panic!("not supported")
        }
   }
}

impl<T: Integer + Montgomery + Clone + Shr<u32, Output = T>> Shr<u32> for Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone, {
    type Output = Self;
    fn shr(self, rhs: u32) -> Self::Output {
        match self.0 {
            Left(v) => Self(Left(v >> rhs)),
            Right(_) => panic!("not supported")
        }
    }
}
impl<T: Integer + Montgomery + Clone + Shr<u32, Output = T>> Shr<u32> for &Mint<T>
where
    T::Double: From<T>,
    T::Inv: Clone, {
    type Output = Mint<T>;
    fn shr(self, rhs: u32) -> Self::Output {
        match &self.0 {
            Left(v) => Mint(Left(v.clone() >> rhs)),
            Right(_) => panic!("not supported")
        }
    }
}



// TODO: implement ModularRefOps

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basics() {
        let a: Mint<u32> = 12.into();
        let b: Mint<u32> = 8.into();
        assert_eq!(a + b, 20.into());
        
        dbg!(crate::nt_funcs::is_prime(&a, None));
    }
}
