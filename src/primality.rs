use std::iter::successors;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{Zero, One, RefNum, NumRef, FromPrimitive};
use std::convert::TryFrom;
use crate::traits::{ModInt, ExactRoots, PrimalityUtils};

/// Utilities for the Lucas pseudoprime test
// TODO (v0.0.5): implement this automatically for ModInt, need a bit iterator trait
pub trait LucasUtils {
    /// Find Lucas sequence to n with modulo m, i.e. $U_{n+1}(P,Q)$ mod m and $V_{n+1}(P,Q)$
    /// Reference: <https://en.wikipedia.org/wiki/Lucas_sequence>
    ///            Peter Hackman, "Elementary Number Theory", section "L.XVII Scalar" <http://hackmat.se/kurser/TATM54/booktot.pdf>
    fn lucasm(p: usize, q: isize, m: Self, n: Self) -> (Self, Self) where Self: Sized;

    /// Find proper parameters P and Q of Lucas sequence for n, such that Jacobi(D|n) is -1
    /// using Selfridge's method A.
    fn pq_selfridge(n: &Self) -> (usize, isize);
}

macro_rules! impl_lucas_util_for_uprim {
    ($T:ty) => {
        impl LucasUtils for $T {
            fn lucasm(p: usize, q: isize, m: Self, n: Self) -> (Self, Self) {
                let p = <$T>::try_from(p).unwrap() % m;
                let q = if q >= 0 { <$T>::try_from(q).unwrap() % &m }
                        // FIXME (v0.1): implement ModInt for <T,T> and prevent this ugly ModInt::<T,T>::xxx call
                        else { ModInt::<&$T, &$T>::negm(&<$T>::try_from(-q).unwrap(), &m) };

                let mut uk = 0; // U(k)
                let mut uk1 = 1; // U(k+1)

                for i in (0..<$T>::BITS).rev() {
                    if (n & (1 << i)) > 0 {
                        // k' = 2k+1
                        // U(k'+1) = U(2k+2) = PU(k+1)² - 2*QU(k+1)U(k)
                        let t1 = p.mulm(&uk1,&m).mulm(&uk1,&m);
                        let t2 = 2.mulm(&q,&m).mulm(&uk1,&m).mulm(&uk,&m);
                        let new_uk1 = t1.subm(t2,&m);
                        // U(k') = U(2k+1) = U(k+1)² - QU(k)²
                        let t1 = uk1.mulm(&uk1,&m);
                        let t2 = q.mulm(&uk,&m).mulm(uk,&m);
                        let new_uk = t1.subm(t2,&m);
                        uk1 = new_uk1; uk = new_uk;
                    } else {
                        // k' = 2k
                        // U(k'+1) = U(2k+1) = U(k+1)² - QU(k)²
                        let t1 = uk1.mulm(&uk1,&m);
                        let t2 = q.mulm(&uk,&m).mulm(&uk,&m);
                        let new_uk1 = t1.subm(t2,&m);
                        // U(k') = U(2k) = 2U(k+1)U(k) - PU(k)²
                        let t1 = 2.mulm(uk1,&m).mulm(&uk,&m);
                        let t2 = p.mulm(&uk,&m).mulm(uk,&m);
                        let new_uk = t1.subm(t2,&m);
                        uk1 = new_uk1; uk = new_uk;
                    }
                }

                let vk = 2.mulm(uk1,&m).subm(p.mulm(&uk,&m),&m);
                (uk, vk)
            }

            fn pq_selfridge(n: &$T) -> (usize, isize) {
                for (d, neg) in successors(Some((5 as $T, false)), |(d, neg)| Some((d+2, !neg))) {
                    let sd = if neg { ModInt::<&$T, &$T>::negm(&d, n) } else { d };
                    if d == 13 && (*n).is_square() { // check if n is a square number after several trial
                        return (0, 0);
                    }
                    if ModInt::<&$T, &$T>::jacobi(&sd, n) == -1 {
                        let d = if neg { -(d as isize) } else { d as isize };
                        return (1, (1 - d as isize) / 4);
                    }
                }
                panic!(); // impossible
            }
        }
    };
}

impl_lucas_util_for_uprim!(u8);
impl_lucas_util_for_uprim!(u16);
impl_lucas_util_for_uprim!(u32);
impl_lucas_util_for_uprim!(u64);
impl_lucas_util_for_uprim!(u128);

impl LucasUtils for BigUint {
    fn lucasm(p: usize, q: isize, m: Self, n: Self) -> (Self, Self) {
        let p = BigUint::from_usize(p).unwrap() % &m;
        let q = if q >= 0 { BigUint::from_isize(q).unwrap() % &m }
                        else { ModInt::<&BigUint, &BigUint>::negm(&BigUint::from_isize(-q).unwrap(), &m) };

        let mut uk = BigUint::zero(); // U(k)
        let mut uk1 = BigUint::one(); // U(k+1)
        let two = BigUint::one() + BigUint::one();

        for i in (0..n.bits()).rev() {
            if n.bit(i) {
                // k' = 2k+1
                // U(k'+1) = U(2k+2) = PU(k+1)² - 2*QU(k+1)U(k)
                let t1 = p.mulm(&uk1,&m).mulm(&uk1,&m);
                let t2 = two.mulm(&q,&m).mulm(&uk1,&m).mulm(&uk,&m);
                let new_uk1 = t1.subm(t2,&m);
                // U(k') = U(2k+1) = U(k+1)² - QU(k)²
                let t1 = uk1.mulm(&uk1,&m);
                let t2 = q.mulm(&uk,&m).mulm(uk,&m);
                let new_uk = t1.subm(t2,&m);
                uk1 = new_uk1; uk = new_uk;
            } else {
                // k' = 2k
                // U(k'+1) = U(2k+1) = U(k+1)² - QU(k)²
                let t1 = uk1.mulm(&uk1,&m);
                let t2 = q.mulm(&uk,&m).mulm(&uk,&m);
                let new_uk1 = t1.subm(t2,&m);
                // U(k') = U(2k) = 2U(k+1)U(k) - PU(k)²
                let t1 = two.mulm(uk1,&m).mulm(&uk,&m);
                let t2 = p.mulm(&uk,&m).mulm(uk,&m);
                let new_uk = t1.subm(t2,&m);
                uk1 = new_uk1; uk = new_uk;
            }
        }

        let vk = two.mulm(uk1,&m).subm(p.mulm(&uk,&m),&m);
        (uk, vk)
    }
    
    fn pq_selfridge(n: &BigUint) -> (usize, isize) {
        unimplemented!();
    }
}

impl<T: Integer + NumRef + Clone + FromPrimitive + LucasUtils> PrimalityUtils for T
where for<'r> &'r T: RefNum<T> + std::ops::Shr<usize, Output = T> + ModInt<&'r T, &'r T, Output = T>
{
    fn is_prp(&self, base: Self) -> bool {
        if self < &Self::one() { return false; }
        let tm1 = self - Self::one();
        base.powm(&tm1, self).is_one()
    }

    fn is_sprp(&self, base: T) -> bool {
        if self < &Self::one() { return false; }

        // find 2^shift*u + 1 = n
        let tm1 = self - T::one();
        let shift = tm1.fac2();
        let u = &tm1 >> shift;

        let mut x = base.powm(&u, self);
        if x == T::one() || x == tm1 { return true }

        for _ in 0..shift {
            x = (&x).mulm(&x, self);
            if x == tm1 { return true }
        }

        x == T::one()
    }

    fn is_lprp(&self, p: Option<usize>, q: Option<isize>) -> bool {
        if self < &Self::one() { return false; }
        if self.is_even() { return false; }

        let (p, q) = match (p, q) {
            (Some(sp), Some(sq)) => (sp, sq),
            (_, _) =>  {
                let (sp, sq) = LucasUtils::pq_selfridge(self);
                if sp == 0 { return false }; // is a perfect power
                (sp, sq)
            }
        };

        let d = (p*p) as isize - 4*q;
        let d = if d > 0 { Self::from_isize(d).unwrap() }
                   else { ModInt::<&Self, &Self>::negm(&Self::from_isize(-d).unwrap(), self) };
        let delta = match d.jacobi(self) {
            0 => self.clone(), -1 => self + Self::one(), 1 => self - Self::one(), _ => panic!("")
        };
        let (u, _) = LucasUtils::lucasm(p, q, self.clone(), delta);
        u.is_zero()

        // TODO: add additional check V(n+1) == 2Q mod n, if p,q are not specified and jacobi = -1
    }

    fn is_slprp(&self, p: Option<usize>, q: Option<isize>) -> bool {
        if self < &Self::one() { return false; }
        if self.is_even() { return false; }

        let (p, q) = match (p, q) {
            (Some(sp), Some(sq)) => (sp, sq),
            (_, _) =>  {
                let (sp, sq) = LucasUtils::pq_selfridge(self);
                if sp == 0 { return false }; // is a perfect power
                (sp, sq)
            }
        };

        let d = (p*p) as isize - 4*q;
        let d = if d > 0 { Self::from_isize(d).unwrap() }
                   else { ModInt::<&Self, &Self>::negm(&Self::from_isize(-d).unwrap(), self) };
        let delta = match d.jacobi(self) {
            0 => self.clone(), -1 => self + Self::one(), 1 => self - Self::one(), _ => panic!("")
        };

        let shift = delta.fac2();
        let base = &delta >> shift;

        let (ud, vd) = LucasUtils::lucasm(p, q, self.clone(), base.clone());
        if ud.is_zero() || vd.is_zero() { return true; }

        // do first iteration to remove the sign on Q
        if shift == 0 { return false; }
        let qk = Self::from_isize(q.abs()).unwrap().powm(&base, self);
        let two = Self::from_u8(2).unwrap();
        // V(2k) = V(k)² - 2Q^k
        let mut vd = if q >= 0 { vd.mulm(&vd,self).subm(&two.mulm(&qk,self),self) }
                        else { vd.mulm(&vd,self).addm(&two.mulm(&qk,self),self) };
        if vd.is_zero() { return true }
        let mut qk = qk.mulm(&qk,self);

        for _ in 1..shift {
            // V(2k) = V(k)² - 2Q^k
            vd = vd.mulm(&vd,self).subm(&two.mulm(&qk,self),self);
            if vd.is_zero() { return true }
            qk = qk.mulm(&qk,self);
        }
        false
    }

    fn is_eslprp(&self, p: Option<usize>) -> bool {
        if self < &Self::one() { return false; }
        if self.is_even() { return false; }

        let p = p.unwrap(); // TODO: find p=3,4,5 such that jacobi = -1, check square if p == 15

        let d = (p*p) as isize - 4;
        let d = if d > 0 { Self::from_isize(d).unwrap() }
                   else { ModInt::<&Self, &Self>::negm(&Self::from_isize(-d).unwrap(), self) };
        let delta = match d.jacobi(self) {
            0 => self.clone(), -1 => self + Self::one(), 1 => self - Self::one(), _ => panic!("")
        };

        let shift = delta.fac2();
        let base = &delta >> shift;

        let (ud, mut vd) = LucasUtils::lucasm(p, 1, self.clone(), base.clone());
        let two = Self::from_u8(2).unwrap();
        // U(d) = 0 or V(d) = ±2
        if ud.is_zero() && (vd == two || vd == self - &two) { return true; }
        if vd.is_zero() { return true; }

        for _ in 1..(shift-1) {
            // V(2k) = V(k)² - 2
            vd = vd.mulm(&vd, self).subm(&two, self);
            if vd.is_zero() { return true }
        }
        false
    }
}

// TODO: Add test using https://github.com/dignifiedquire/num-bigint/blob/master/src/prime.rs

#[cfg(test)]
mod tests {
    use super::*;
    use rand::random;

    #[test]
    fn lucas_mod_test() {
        // OEIS A006190
        let p3qm1: [u64; 26] = [0, 1, 3, 10, 33, 109, 360, 1189, 3927, 12970, 42837, 141481, 467280, 1543321, 5097243, 16835050, 55602393, 183642229, 606529080, 2003229469, 6616217487, 21851881930, 72171863277, 238367471761, 787274278560, 2600190307441];
        let m = random::<u16>();
        for n in 2..p3qm1.len() {
            let (uk, _) = LucasUtils::lucasm(3, -1, m as u64, n as u64);
            assert_eq!(uk, p3qm1[n] % (m as u64));
            let (uk, _) = LucasUtils::lucasm(3, -1, BigUint::from(m), BigUint::from(n));
            assert_eq!(uk, BigUint::from(p3qm1[n] % (m as u64)));
        }
    }

    // TODO: add test for is_prp and is_sprp

    #[test]
    fn lucas_prp_test() {
        assert!(19u8.is_lprp(Some(3), Some(-1)));

        // least lucas pseudo primes for Q=-1 and Jacobi(D/n) = -1 (from Wikipedia)
        let plimit: [u16; 5] = [323, 35, 119, 9, 9];
        for (i, l) in plimit.iter().cloned().enumerate() {
            let p = i + 1;
            assert!(l.is_lprp(Some(p), Some(-1)));

            // test four random numbers under limit
            for _ in 0..10 {
                let n = random::<u16>() % l;
                if n < 2 || n.is_sprp(2) { continue } // skip real primes
                let d = (p*p + 4) as u16;
                if n.is_odd() && ModInt::<&u16,&u16>::jacobi(&d, &n) != -1 { continue }
                assert!(!n.is_lprp(Some(p), Some(-1)), "lucas prp test on {} with p = {}", n, p);
            }
        }

        // least strong lucas pseudoprimes for Q=-1 and Jacobi(D/n) = -1 (from Wikipedia)
        let plimit: [u16; 3] = [4181, 169, 119];
        for (i, l) in plimit.iter().cloned().enumerate() {
            let p = i + 1;
            assert!(l.is_slprp(Some(p), Some(-1)));

            // test random numbers under limit
            for _ in 0..10 {
                let n = random::<u16>() % l;
                if n < 2 || (n.is_sprp(2) && n.is_sprp(3)) { continue } // skip real primes
                let d = (p*p + 4) as u16;
                if n.is_odd() && ModInt::<&u16,&u16>::jacobi(&d, &n) != -1 { continue }
                assert!(!n.is_slprp(Some(p), Some(-1)), "strong lucas prp test on {} with p = {}", n, p);
            }
        }

        // lucas pseudoprimes under 10000, OEIS A217120
        let lpsp: [u16; 9] = [323, 377, 1159, 1829, 3827, 5459, 5777, 9071, 9179];
        for psp in lpsp {
            assert!(psp.is_lprp(None, None), "lucas prp test on {}", psp);
        }
        for _ in 0..10 {
            let n = random::<u16>() % 10000;
            if n < 2 || (n.is_sprp(2) && n.is_sprp(3)) { continue } // skip real primes
            if lpsp.iter().find(|&x| x == &n).is_some() { continue } // skip pseudo primes
            assert!(!n.is_lprp(None, None), "lucas prp test on {}", n);
        }
        
        // strong lucas pseudoprimes under 10000, OEIS A217255
        let slpsp: [u16; 2] = [5459, 5777];
        for psp in slpsp {
            assert!(psp.is_slprp(None, None), "strong lucas prp test on {}", psp);
        }
        for _ in 0..10 {
            let n = random::<u16>() % 10000;
            if n < 2 || (n.is_sprp(2) && n.is_sprp(3)) { continue } // skip real primes
            if slpsp.iter().find(|&x| x == &n).is_some() { continue } // skip pseudo primes
            assert!(!n.is_slprp(None, None), "strong lucas prp test on {}", n);
        }
    }
}
