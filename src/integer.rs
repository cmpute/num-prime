// backend implementations for integers
use num_bigint::BigUint;
use num_traits::{ToPrimitive};
use crate::traits::{ModInt};

macro_rules! impl_mod_arithm_uu {
    ($T:ty, $Tdouble:ty) => {
        impl ModInt<&$T, &$T> for &$T {
            type Output = $T;
            #[inline]
            fn fac2(self) -> usize { self.trailing_zeros() as usize }
            fn mulm(self, rhs: &$T, m: &$T) -> $T { ((*self as $Tdouble) * (*rhs as $Tdouble) % (*m as $Tdouble)) as $T }
            fn powm(self, exp: &$T, m: &$T) -> $T {
                if *exp == 1 { return self % m; }
                if *exp == 2 { return self.mulm(self, m); }
        
                let mut multi = self % m;
                let mut exp = *exp;
                let mut result = 1;
                while exp > 0 {
                    if exp & 1 > 0 {
                        result = result.mulm(&multi, m);
                    }
                    multi = multi.mulm(&multi, m);
                    exp >>= 1;
                }
                result
            }
        }
        
        impl ModInt<$T, &$T> for &$T {
            type Output = $T;
            #[inline]
            fn fac2(self) -> usize { self.trailing_zeros() as usize }
            #[inline]
            fn mulm(self, rhs: $T, m: &$T) -> $T { self.mulm(&rhs, m) }
            #[inline]
            fn powm(self, exp: $T, m: &$T) -> $T { self.powm(&exp, m) }
        }
    }
}

impl_mod_arithm_uu!(u8, u16);
impl_mod_arithm_uu!(u16, u32);
impl_mod_arithm_uu!(u32, u64);
impl_mod_arithm_uu!(u64, u128);

impl ModInt<&u128, &u128> for &u128 {
    type Output = u128;

    #[inline]
    fn fac2(self) -> usize { self.trailing_zeros() as usize }

    // TODO: check if this is also faster in u64
    // TODO: benchmark against http://www.janfeitsma.nl/math/psp2/expmod
    fn mulm(self, rhs: &u128, m: &u128) -> u128 {
        if let Some(ab) = self.checked_mul(*rhs) {
            return ab % m
        }

        let mut a = self % m;
        let mut b = rhs % m;

        if let Some(ab) = a.checked_mul(b) {
            return ab % m
        }

        let mut result: u128 = 0;
        while b > 0 {
            if b & 1 > 0 {
                result += a;
                result %= m;
            }
            a <<= 1;
            if a >= *m {
                a %= m;
            }
            b >>= 1;
        }
        result
    }

    fn powm(self, exp: &u128, m: &u128) -> u128 {
        if *exp == 1 {
            return self % m;
        }

        let mut multi = self % m;
        let mut exp = *exp;
        let mut result = 1;
        while exp > 0 {
            if exp & 1 > 0 {
                result = result.mulm(&multi, m);
            }
            multi = multi.mulm(&multi, m);
            exp >>= 1;
        }
        result
    }
}

impl ModInt<&BigUint, &BigUint> for &BigUint {    
    type Output = BigUint;

    #[inline]
    fn fac2(self) -> usize { 
        match BigUint::trailing_zeros(self) {
            Some(a) => a as usize, None => 0
        }
    }

    fn mulm(self, rhs: &BigUint, m: &BigUint) -> BigUint {
        let a = self % m;
        let b = rhs % m;

        if let Some(sm) = m.to_u64() {
            let sself = a.to_u64().unwrap();
            let srhs = b.to_u64().unwrap();
            return BigUint::from(sself.mulm(&srhs, &sm));
        }

        (a * b) % m
    }

    #[inline]
    fn powm(self, exp: &BigUint, m: &BigUint) -> BigUint {
        self.modpow(&exp, m)
    }
}

impl ModInt<BigUint, &BigUint> for &BigUint {
    type Output = BigUint;
    
    #[inline]
    fn fac2(self) -> usize { 
        match BigUint::trailing_zeros(self) {
            Some(a) => a as usize, None => 0
        }
    }
    #[inline]
    fn mulm(self, rhs: BigUint, m: &BigUint) -> BigUint { self.mulm(&rhs, m) }
    #[inline]
    fn powm(self, exp: BigUint, m: &BigUint) -> BigUint { self.powm(&exp, m) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand;
    use num_bigint::RandBigInt;

    #[test]
    fn u64_mod_test() {
        let a = rand::random::<u64>() % 100000;
        let m = rand::random::<u64>() % 100000;
        assert_eq!(a.mulm(a, &m), (a * a) % m);
        assert_eq!(a.powm(3, &m), a.pow(3) % m);
    }

    #[test]
    fn biguint_mod_test() {
        let mut rng = rand::thread_rng();
        let a = rng.gen_biguint(500); let ra = &a;
        let m = rng.gen_biguint(500); let rm = &m;
        assert_eq!(ra.mulm(ra, rm), (ra * ra) % rm);
        assert_eq!(ra.powm(BigUint::from(3u8), rm), ra.pow(3) % rm);
    }
}
