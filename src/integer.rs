//! Backend implementations for integers

use crate::traits::{BitTest, ExactRoots};
use num_integer::Integer;

#[cfg(feature = "num-bigint")]
use num_bigint::BigUint;
#[cfg(feature = "num-bigint")]
use num_traits::{One, ToPrimitive, Zero};

macro_rules! impl_bititer_prim {
    ($($T:ty)*) => {$(
        impl BitTest for $T {
            #[inline]
            fn bits(&self) -> usize {
                (<$T>::BITS - self.leading_zeros()) as usize
            }
            #[inline]
            fn bit(&self, position: usize) -> bool {
                self & (1 << position) > 0
            }
            #[inline]
            fn trailing_zeros(&self) -> usize {
                <$T>::trailing_zeros(*self) as usize
            }
        }
    )*}
}
impl_bititer_prim!(u8 u16 u32 u64 u128 usize);

#[cfg(feature = "num-bigint")]
impl BitTest for BigUint {
    fn bit(&self, position: usize) -> bool {
        self.bit(position as u64)
    }
    fn bits(&self) -> usize {
        BigUint::bits(&self) as usize
    }
    #[inline]
    fn trailing_zeros(&self) -> usize {
        match BigUint::trailing_zeros(&self) {
            Some(a) => a as usize,
            None => 0,
        }
    }
}

// QUAD_RESIDUAL[N] has a bit i set iff i is a quadratic residue mod N.
const QUAD_RESIDUAL64: u64 = 0x0202021202030213;
const QUAD_RESIDUAL63: u64 = 0x0402483012450293;
const QUAD_RESIDUAL65: u64 = 0x218a019866014613;
const QUAD_RESIDUAL11: u64 = 0x23b;

macro_rules! impl_exactroot_prim {
    ($($T:ty)*) => {$(
        impl ExactRoots for $T {
            fn sqrt_exact(&self) -> Option<Self> {
                // eliminate most non-squares by checking legendre symbols.
                // See H. Cohen's "Course in Computational Algebraic Number Theory",
                // algorithm 1.7.3, page 40.
                if (QUAD_RESIDUAL64 >> (self & 63)) & 1 == 0 {
                    return None;
                }
                if (QUAD_RESIDUAL63 >> (self % 63)) & 1 == 0 {
                    return None;
                }
                if (QUAD_RESIDUAL65 >> ((self % 65) & 63)) & 1 == 0 {
                    // Both 0 and 64 are squares mod 65
                    return None;
                }
                if (QUAD_RESIDUAL11 >> (self % 11)) & 1 == 0 {
                    return None;
                }
                self.nth_root_exact(2)
            }
        }
    )*};
}
impl_exactroot_prim!(u8 u16 u32 u64 u128 usize);

#[cfg(feature = "num-bigint")]
impl ExactRoots for BigUint {
    // TODO: improve based on SqrtExact.java @ java-math-library
    fn sqrt_exact(&self) -> Option<Self> {
        if (QUAD_RESIDUAL64 >> (self % 64u8).to_u64().unwrap()) & 1 == 0 {
            return None;
        }
        if (QUAD_RESIDUAL63 >> (self % 63u8).to_u64().unwrap()) & 1 == 0 {
            return None;
        }
        if (QUAD_RESIDUAL65 >> ((self % 65u8) % 64u8).to_u64().unwrap()) & 1 == 0 {
            return None;
        }
        if (QUAD_RESIDUAL11 >> (self % 11u8).to_u64().unwrap()) & 1 == 0 {
            return None;
        }
        self.nth_root_exact(2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand;

    #[test]
    fn exact_root_test() {
        // some simple tests
        assert!(matches!(ExactRoots::sqrt_exact(&3u8), None));
        assert!(matches!(ExactRoots::sqrt_exact(&4u8), Some(2)));
        assert!(matches!(ExactRoots::sqrt_exact(&9u8), Some(3)));
        assert!(matches!(ExactRoots::sqrt_exact(&18u8), None));

        // test fast implementations of sqrt against nth_root
        for _ in 0..100 {
            let x = rand::random::<u32>();
            assert_eq!(
                ExactRoots::sqrt_exact(&x),
                ExactRoots::nth_root_exact(&x, 2)
            );
        }
        for _ in 0..100 {
            let x = rand::random::<u32>() as u64;
            assert!(matches!(ExactRoots::sqrt_exact(&(x * x)), Some(v) if v == x));
        }
        for _ in 0..100 {
            let x = rand::random::<u32>() as u64;
            let y = rand::random::<u32>() as u64;
            if x == y {
                continue;
            }
            assert!(ExactRoots::sqrt_exact(&(x * y)).is_none());
        }
    }
}
