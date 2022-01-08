// backend implementations for integers
use num_bigint::BigUint;
use num_traits::{ToPrimitive};
use crate::traits::{ModInt, ArithmeticHelpers};

impl ArithmeticHelpers for u64 {
    /// Returns greatest common divisor between a, b
    #[inline]
    fn trailing_zeros(&self) -> usize { u64::trailing_zeros(*self) as usize }
}

impl ModInt<&u64, &u64> for &u64 {
    type Output = u64;

    fn mul_mod(self, rhs: &u64, m: &u64) -> u64 {
        if let Some(ab) = self.checked_mul(*rhs) {
            return ab % m
        }

        let mut a = self % m;
        let mut b = rhs % m;

        if let Some(ab) = a.checked_mul(b) {
            return ab % m
        }

        let mut result: u64 = 0;
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

    fn pow_mod(self, exp: &u64, m: &u64) -> u64 {
        if *exp == 1 {
            return self % m;
        }

        if *exp < (u32::MAX as u64) {
            if let Some(ae) = self.checked_pow(*exp as u32) {
                return ae % m;
            }
        }

        let mut multi = self % m;
        let mut exp = *exp;
        let mut result = 1;
        while exp > 0 {
            if exp & 1 > 0 {
                result = result.mul_mod(&multi, m);
            }
            multi = multi.mul_mod(&multi, m);
            exp >>= 1;
        }
        result
    }
}

impl ModInt<u64, &u64> for &u64 {
    type Output = u64;
    #[inline]
    fn mul_mod(self, rhs: u64, m: &u64) -> u64 { self.mul_mod(&rhs, m) }
    #[inline]
    fn pow_mod(self, exp: u64, m: &u64) -> u64 { self.pow_mod(&exp, m) }
}

impl ArithmeticHelpers for BigUint {
    /// Returns greatest common divisor between a, b
    #[inline]
    fn trailing_zeros(&self) -> usize { 
        match BigUint::trailing_zeros(&self) {
            Some(a) => a as usize, None => 0
        }
    }
}

impl ModInt<&BigUint, &BigUint> for &BigUint {    
    type Output = BigUint;

    fn mul_mod(self, rhs: &BigUint, m: &BigUint) -> BigUint {
        let a = self % m;
        let b = rhs % m;

        if let Some(sm) = m.to_u64() {
            let sself = a.to_u64().unwrap();
            let srhs = b.to_u64().unwrap();
            return BigUint::from(sself.mul_mod(&srhs, &sm));
        }

        (a * b) % m
    }

    #[inline]
    fn pow_mod(self, exp: &BigUint, m: &BigUint) -> BigUint {
        self.modpow(&exp, m)
    }
}

impl ModInt<BigUint, &BigUint> for BigUint {
    type Output = BigUint;
    #[inline]
    fn mul_mod(self, rhs: BigUint, m: &BigUint) -> BigUint { self.mul_mod(&rhs, m) }
    #[inline]
    fn pow_mod(self, exp: BigUint, m: &BigUint) -> BigUint { self.pow_mod(&exp, m) }
}
impl ModInt<BigUint, &BigUint> for &BigUint {
    type Output = BigUint;
    #[inline]
    fn mul_mod(self, rhs: BigUint, m: &BigUint) -> BigUint { self.mul_mod(&rhs, m) }
    #[inline]
    fn pow_mod(self, exp: BigUint, m: &BigUint) -> BigUint { self.pow_mod(&exp, m) }
}
impl ModInt<&BigUint, &BigUint> for BigUint {
    type Output = BigUint;
    #[inline]
    fn mul_mod(self, rhs: &BigUint, m: &BigUint) -> BigUint { (&self).mul_mod(rhs, m) }
    #[inline]
    fn pow_mod(self, exp: &BigUint, m: &BigUint) -> BigUint { (&self).pow_mod(exp, m) }
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
        assert_eq!(a.mul_mod(a, &m), (a * a) % m);
        assert_eq!(a.pow_mod(3, &m), a.pow(3) % m);
    }

    #[test]
    fn biguint_mod_test() {
        let mut rng = rand::thread_rng();
        let a = rng.gen_biguint(500); let ra = &a;
        let m = rng.gen_biguint(500); let rm = &m;
        assert_eq!(ra.mul_mod(ra, rm), (ra * ra) % rm);
        assert_eq!(ra.pow_mod(BigUint::from(3u8), rm), ra.pow(3) % rm);
    }
}
