//! Backend implementations for integers

use crate::traits::{ExactRoots, ModInt, BitTest};
use num_bigint::BigUint;
use num_integer::{Integer, Roots};
use num_traits::{One, Pow, ToPrimitive, Zero};


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

// TODO (v0.1): implement fast perfect power check for integers
// REF: https://github.com/coreutils/coreutils/blob/master/src/factor.c#L1833
//      https://math.stackexchange.com/a/878338
//      https://github.com/nbraud/coreutils/blob/factor/sqrt/src/uu/factor/src/numeric/sqrt.rs#L35
impl<T: Roots + Pow<u32, Output = Self> + Clone> ExactRoots for T {}

macro_rules! impl_jacobi_prim {
    ($T:ty) => {
        fn jacobi(self, n: &$T) -> i8 {
            if n % 2 == 0 || n < &0 {
                panic!("The Jacobi symbol is only defined for non-negative odd integers!")
            }

            if self == 0 {
                return 0;
            }
            if self == 1 {
                return 1;
            }

            let mut a = self % n;
            let mut n = n.clone();
            let mut t = 1;
            while a > 0 {
                while (a & 1) == 0 {
                    a = a / 2;
                    if n & 7 == 3 || n & 7 == 5 {
                        t *= -1;
                    }
                }
                std::mem::swap(&mut a, &mut n);
                if (a & 3) == 3 && (n & 3) == 3 {
                    t *= -1;
                }
                a = a % n;
            }
            if n == 1 {
                t
            } else {
                0
            }
        }
    };
}

// implement inverse mod using extended euclidean algorithm
macro_rules! impl_invm_prim {
    ($T:ty) => {
        fn invm(self, m: &$T) -> Option<Self::Output> {
            let x = if &self >= m { self % m } else { self.clone() };

            let (mut last_r, mut r) = (m.clone(), x);
            let (mut last_t, mut t) = (0, 1);

            while r > 0 {
                let (quo, rem) = last_r.div_rem(&r);
                last_r = r;
                r = rem;

                let new_t = last_t.subm(quo.mulm(t, m), m);
                last_t = t;
                t = new_t;
            }

            // if r = gcd(self, m) > 1, then inverse doesn't exist
            if last_r > 1 {
                None
            } else {
                Some(last_t)
            }
        }
    };
}

macro_rules! impl_mod_arithm_uu {
    ($T:ty, $Tdouble:ty) => {
        impl ModInt<$T, &$T> for $T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: $T, m: &$T) -> $T {
                (((self as $Tdouble) + (rhs as $Tdouble)) % (*m as $Tdouble)) as $T
            }
            #[inline]
            fn subm(self, rhs: $T, m: &$T) -> $T {
                let (lhs, rhs) = (self % m, rhs % m);
                if lhs >= rhs {
                    lhs - rhs
                } else {
                    m - (rhs - lhs)
                }
            }
            #[inline]
            fn mulm(self, rhs: $T, m: &$T) -> $T {
                (((self as $Tdouble) * (rhs as $Tdouble)) % (*m as $Tdouble)) as $T
            }
            fn powm(self, exp: $T, m: &$T) -> $T {
                if exp == 1 {
                    return self % m;
                }
                if exp == 2 {
                    return self.mulm(self, m);
                }

                let mut multi = self % m;
                let mut exp = exp;
                let mut result = 1;
                while exp > 0 {
                    if exp & 1 > 0 {
                        result = result.mulm(multi, m);
                    }
                    multi = multi.mulm(multi, m);
                    exp >>= 1;
                }
                result
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                let x = self % m;
                if x == 0 {
                    0
                } else {
                    m - x
                }
            }
            impl_jacobi_prim!($T);
            impl_invm_prim!($T);
        }
    };
}

impl_mod_arithm_uu!(u8, u16);
impl_mod_arithm_uu!(u16, u32);
impl_mod_arithm_uu!(u32, u64);
impl_mod_arithm_uu!(u64, u128);

impl ModInt<u128, &u128> for u128 {
    type Output = u128;

    // XXX: check if these operations are also faster in u64
    #[inline]
    fn addm(self, rhs: u128, m: &u128) -> u128 {
        if let Some(ab) = self.checked_add(rhs) {
            return ab % m;
        }

        let (lhs, rhs) = (self % m, rhs % m);
        if lhs < m - rhs {
            lhs + rhs
        } else {
            lhs.min(rhs) - (m - lhs.max(rhs))
        }
    }

    #[inline]
    fn subm(self, rhs: u128, m: &u128) -> u128 {
        let (lhs, rhs) = (self % m, rhs % m);
        if lhs >= rhs {
            lhs - rhs
        } else {
            m - (rhs - lhs)
        }
    }

    // TODO: benchmark against http://www.janfeitsma.nl/math/psp2/expmod
    fn mulm(self, rhs: u128, m: &u128) -> u128 {
        if let Some(ab) = self.checked_mul(rhs) {
            return ab % m;
        }

        let mut a = self % m;
        let mut b = rhs % m;

        if let Some(ab) = a.checked_mul(b) {
            return ab % m;
        }

        let mut result: u128 = 0;
        while b > 0 {
            if b & 1 > 0 {
                result = result.addm(a, m);
            }
            a = a.addm(a, m);
            b >>= 1;
        }
        result
    }

    fn powm(self, exp: u128, m: &u128) -> u128 {
        if exp == 1 {
            return self % m;
        }

        let mut multi = self % m;
        let mut exp = exp;
        let mut result = 1;
        while exp > 0 {
            if exp & 1 > 0 {
                result = result.mulm(multi, m);
            }
            multi = multi.mulm(multi, m);
            exp >>= 1;
        }
        result
    }

    fn negm(self, m: &u128) -> u128 {
        let x = self % m;
        if x == 0 {
            0
        } else {
            m - x
        }
    }

    impl_jacobi_prim!(u128);
    impl_invm_prim!(u128);
}

macro_rules! impl_mod_arithm_by_deref {
    ($T:ty) => {
        impl ModInt<$T, &$T> for &$T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: $T, m: &$T) -> $T {
                (*self).addm(rhs, &m)
            }
            #[inline]
            fn subm(self, rhs: $T, m: &$T) -> $T {
                (*self).subm(rhs, &m)
            }
            #[inline]
            fn mulm(self, rhs: $T, m: &$T) -> $T {
                (*self).mulm(rhs, &m)
            }
            #[inline]
            fn powm(self, exp: $T, m: &$T) -> $T {
                (*self).powm(exp, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModInt::<$T, &$T>::negm(*self, m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModInt::<$T, &$T>::invm(*self, m)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModInt::<$T, &$T>::jacobi(*self, n)
            }
        }

        impl ModInt<&$T, &$T> for $T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: &$T, m: &$T) -> $T {
                self.addm(*rhs, &m)
            }
            #[inline]
            fn subm(self, rhs: &$T, m: &$T) -> $T {
                self.subm(*rhs, &m)
            }
            #[inline]
            fn mulm(self, rhs: &$T, m: &$T) -> $T {
                self.mulm(*rhs, &m)
            }
            #[inline]
            fn powm(self, exp: &$T, m: &$T) -> $T {
                self.powm(*exp, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModInt::<$T, &$T>::negm(self, m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModInt::<$T, &$T>::invm(self, m)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModInt::<$T, &$T>::jacobi(self, n)
            }
        }

        impl ModInt<&$T, &$T> for &$T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: &$T, m: &$T) -> $T {
                (*self).addm(*rhs, &m)
            }
            #[inline]
            fn subm(self, rhs: &$T, m: &$T) -> $T {
                (*self).subm(*rhs, &m)
            }
            #[inline]
            fn mulm(self, rhs: &$T, m: &$T) -> $T {
                (*self).mulm(*rhs, &m)
            }
            #[inline]
            fn powm(self, exp: &$T, m: &$T) -> $T {
                (*self).powm(*exp, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModInt::<$T, &$T>::negm(*self, m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModInt::<$T, &$T>::invm(*self, m)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModInt::<$T, &$T>::jacobi(*self, n)
            }
        }
    }
}

impl_mod_arithm_by_deref!(u8);
impl_mod_arithm_by_deref!(u16);
impl_mod_arithm_by_deref!(u32);
impl_mod_arithm_by_deref!(u64);
impl_mod_arithm_by_deref!(u128);

impl ModInt<&BigUint, &BigUint> for &BigUint {
    type Output = BigUint;

    #[inline]
    fn addm(self, rhs: &BigUint, m: &BigUint) -> BigUint {
        (self + rhs) % m
    }
    fn subm(self, rhs: &BigUint, m: &BigUint) -> BigUint {
        let (lhs, rhs) = (self % m, rhs % m);
        if lhs >= rhs {
            lhs - rhs
        } else {
            m - (rhs - lhs)
        }
    }

    fn mulm(self, rhs: &BigUint, m: &BigUint) -> BigUint {
        let a = self % m;
        let b = rhs % m;

        if let Some(sm) = m.to_u64() {
            let sself = a.to_u64().unwrap();
            let srhs = b.to_u64().unwrap();
            return BigUint::from(sself.mulm(srhs, &sm));
        }

        (a * b) % m
    }

    #[inline]
    fn powm(self, exp: &BigUint, m: &BigUint) -> BigUint {
        self.modpow(&exp, m)
    }

    #[inline]
    fn negm(self, m: &BigUint) -> BigUint {
        let x = self % m;
        if x.is_zero() {
            BigUint::zero()
        } else {
            m - x
        }
    }

    fn jacobi(self, n: &BigUint) -> i8 {
        debug_assert!(n.is_odd());

        if self.is_zero() {
            return 0;
        }
        if self.is_one() {
            return 1;
        }

        let three = BigUint::from(3u8);
        let five = BigUint::from(5u8);
        let seven = BigUint::from(7u8);

        let mut a = self % n;
        let mut n = n.clone();
        let mut t = 1;
        while a > BigUint::zero() {
            while a.is_even() {
                a >>= 1;
                if &n & &seven == three || &n & &seven == five {
                    t *= -1;
                }
            }
            std::mem::swap(&mut a, &mut n);
            if (&a & &three) == three && (&n & &three) == three {
                t *= -1;
            }
            a %= &n;
        }
        if n.is_one() {
            t
        } else {
            0
        }
    }

    fn invm(self, m: &BigUint) -> Option<Self::Output> {
        let x = if self >= m { self % m } else { self.clone() };

        let (mut last_r, mut r) = (m.clone(), x);
        let (mut last_t, mut t) = (BigUint::zero(), BigUint::one());

        while r > BigUint::zero() {
            let (quo, rem) = last_r.div_rem(&r);
            last_r = r;
            r = rem;

            let new_t = last_t.subm(&quo.mulm(&t, m), m);
            last_t = t;
            t = new_t;
        }

        // if r = gcd(self, m) > 1, then inverse doesn't exist
        if last_r > BigUint::one() {
            None
        } else {
            Some(last_t)
        }
    }
}

macro_rules! impl_mod_arithm_by_ref {
    ($T:ty) => {
        impl ModInt<$T, &$T> for &$T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: $T, m: &$T) -> $T {
                self.addm(&rhs, &m)
            }
            #[inline]
            fn subm(self, rhs: $T, m: &$T) -> $T {
                self.subm(&rhs, &m)
            }
            #[inline]
            fn mulm(self, rhs: $T, m: &$T) -> $T {
                self.mulm(&rhs, &m)
            }
            #[inline]
            fn powm(self, exp: $T, m: &$T) -> $T {
                self.powm(&exp, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModInt::<&$T, &$T>::negm(self, m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModInt::<&$T, &$T>::invm(self, m)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModInt::<&$T, &$T>::jacobi(self, n)
            }
        }

        impl ModInt<&$T, &$T> for $T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: &$T, m: &$T) -> $T {
                (&self).addm(rhs, &m)
            }
            #[inline]
            fn subm(self, rhs: &$T, m: &$T) -> $T {
                (&self).subm(rhs, &m)
            }
            #[inline]
            fn mulm(self, rhs: &$T, m: &$T) -> $T {
                (&self).mulm(rhs, &m)
            }
            #[inline]
            fn powm(self, exp: &$T, m: &$T) -> $T {
                (&self).powm(exp, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModInt::<&$T, &$T>::negm(&self, m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModInt::<&$T, &$T>::invm(&self, m)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModInt::<&$T, &$T>::jacobi(&self, n)
            }
        }

        impl ModInt<$T, &$T> for $T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: $T, m: &$T) -> $T {
                (&self).addm(&rhs, &m)
            }
            #[inline]
            fn subm(self, rhs: $T, m: &$T) -> $T {
                (&self).subm(&rhs, &m)
            }
            #[inline]
            fn mulm(self, rhs: $T, m: &$T) -> $T {
                (&self).mulm(&rhs, &m)
            }
            #[inline]
            fn powm(self, exp: $T, m: &$T) -> $T {
                (&self).powm(&exp, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModInt::<&$T, &$T>::negm(&self, m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModInt::<&$T, &$T>::invm(&self, m)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModInt::<&$T, &$T>::jacobi(&self, n)
            }
        }
    }
}

impl_mod_arithm_by_ref!(BigUint);

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::RandBigInt;
    use rand;

    #[test]
    fn u64_basic_mod_test() {
        let a = rand::random::<u64>() % 100000;
        let m = rand::random::<u64>() % 100000;
        assert_eq!(a.addm(a, &m), (a + a) % m);
        assert_eq!(a.mulm(a, &m), (a * a) % m);
        assert_eq!(a.powm(3, &m), a.pow(3) % m);
    }

    #[test]
    fn biguint_basic_mod_test() {
        let mut rng = rand::thread_rng();
        let a = rng.gen_biguint(500);
        let ra = &a;
        let m = rng.gen_biguint(500);
        let rm = &m;
        assert_eq!(ra.addm(ra, rm), (ra + ra) % rm);
        assert_eq!(ra.mulm(ra, rm), (ra * ra) % rm);
        assert_eq!(ra.powm(BigUint::from(3u8), rm), ra.pow(3) % rm);
    }

    #[test]
    fn addm_test() {
        let m = 5;

        let test_cases: [(u8, u8, u8); 10] = [
            // [x, y, rem]: x + y = rem (mod m)
            (0, 0, 0),
            (1, 2, 3),
            (2, 1, 3),
            (2, 2, 4),
            (3, 2, 0),
            (2, 3, 0),
            (6, 1, 2),
            (1, 6, 2),
            (11, 7, 3),
            (7, 11, 3),
        ];

        for (x, y, r) in test_cases.iter() {
            assert_eq!(x.addm(y, &m), *r, "u8 x: {}, y: {}", x, y);
            assert_eq!(
                (*x as u16).addm(*y as u16, &(m as u16)),
                *r as u16,
                "u16 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                (*x as u32).addm(*y as u32, &(m as u32)),
                *r as u32,
                "u32 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                (*x as u64).addm(*y as u64, &(m as u64)),
                *r as u64,
                "u64 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                (*x as u128).addm(*y as u128, &(m as u128)),
                *r as u128,
                "u128 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                BigUint::from(*x).addm(BigUint::from(*y), &BigUint::from(m)),
                BigUint::from(*r),
                "biguint x: {}, y: {}",
                x,
                y
            );
        }
    }

    #[test]
    fn subm_test() {
        let m = 7;

        let test_cases: [(u8, u8, u8); 10] = [
            // [x, y, rem]: x - y = rem (mod modu)
            (0, 0, 0),
            (11, 9, 2),
            (5, 2, 3),
            (2, 5, 4),
            (6, 7, 6),
            (1, 7, 1),
            (7, 1, 6),
            (0, 6, 1),
            (15, 1, 0),
            (1, 15, 0),
        ];

        for (x, y, r) in test_cases.iter() {
            assert_eq!(x.subm(y, &m), *r, "u8 x: {}, y: {}", x, y);
            assert_eq!(
                (*x as u16).subm(*y as u16, &(m as u16)),
                *r as u16,
                "u16 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                (*x as u32).subm(*y as u32, &(m as u32)),
                *r as u32,
                "u32 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                (*x as u64).subm(*y as u64, &(m as u64)),
                *r as u64,
                "u64 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                (*x as u128).subm(*y as u128, &(m as u128)),
                *r as u128,
                "u128 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                BigUint::from(*x).subm(BigUint::from(*y), &BigUint::from(m)),
                BigUint::from(*r),
                "biguint x: {}, y: {}",
                x,
                y
            );
        }
    }

    #[test]
    fn invm_test() {
        let test_cases: [(u64, u64, u64); 8] = [
            // [a, m, x] s.t. a*x = 1 (mod m) is satisfied
            (5, 11, 9),
            (8, 11, 7),
            (10, 11, 10),
            (3, 5000, 1667),
            (1667, 5000, 3),
            (999, 5000, 3999),
            (999, 9_223_372_036_854_775_807, 3_619_181_019_466_538_655),
            (
                9_223_372_036_854_775_804,
                9_223_372_036_854_775_807,
                3_074_457_345_618_258_602,
            ),
        ];

        for (a, m, x) in test_cases.iter() {
            assert_eq!(
                ModInt::<&u64>::invm(a, m).unwrap(),
                *x,
                "a: {}, m: {}",
                a,
                m
            );
            assert_eq!(
                ModInt::<&BigUint>::invm(&BigUint::from(*a), &BigUint::from(*m)).unwrap(),
                BigUint::from(*x),
                "a: {}, m: {}",
                a,
                m
            );
        }
    }

    #[test]
    fn jacobi_test() {
        let test_cases: [(u8, u8, i8); 15] = [
            (1, 1, 1),
            (15, 1, 1),
            (2, 3, -1),
            (29, 9, 1),
            (4, 11, 1),
            (17, 11, -1),
            (19, 29, -1),
            (10, 33, -1),
            (11, 33, 0),
            (12, 33, 0),
            (14, 33, -1),
            (15, 33, 0),
            (15, 37, -1),
            (29, 59, 1),
            (30, 59, -1),
        ];

        for (a, n, res) in test_cases.iter() {
            assert_eq!(ModInt::<&u8>::jacobi(a, n), *res, "u8 a: {}, n: {}", a, n);
            assert_eq!(
                ModInt::<&u16>::jacobi(&(*a as u16), &(*n as u16)),
                *res,
                "u16 a: {}, n: {}",
                a,
                n
            );
            assert_eq!(
                ModInt::<&u32>::jacobi(&(*a as u32), &(*n as u32)),
                *res,
                "u32 a: {}, n: {}",
                a,
                n
            );
            assert_eq!(
                ModInt::<&u64>::jacobi(&(*a as u64), &(*n as u64)),
                *res,
                "u64 a: {}, n: {}",
                a,
                n
            );
            assert_eq!(
                ModInt::<&u128>::jacobi(&(*a as u128), &(*n as u128)),
                *res,
                "u128 a: {}, n: {}",
                a,
                n
            );
            assert_eq!(
                ModInt::<&BigUint>::jacobi(&(BigUint::from(*a)), &(BigUint::from(*n))),
                *res,
                "u32 a: {}, n: {}",
                a,
                n
            );
        }
    }
}
