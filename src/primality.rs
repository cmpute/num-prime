use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{Zero, One, RefNum, NumRef, FromPrimitive};
use std::convert::TryFrom;
use crate::traits::{ModInt, PrimalityUtils};

/// Find Lucas sequence to n with modulo m, i.e. $U_{n+1}(P,Q)$ mod m and $V_{n+1}(P,Q)$
/// Reference: <https://en.wikipedia.org/wiki/Lucas_sequence>
///            Peter Hackman, "Elementary Number Theory", section "L.XVII Scalar" <http://hackmat.se/kurser/TATM54/booktot.pdf>
pub trait LucasMod {
    fn lucasm(p: isize, q: isize, m: Self, n: Self) -> (Self, Self) where Self: Sized;
}

macro_rules! impl_lucasm_for_uprim {
    ($T:ty) => {
        impl LucasMod for $T {
            fn lucasm(p: isize, q: isize, m: Self, n: Self) -> (Self, Self) {
                let p = if p >= 0 { <$T>::try_from(p).unwrap() % m }
                                else { ModInt::<&$T, &$T>::negm(&<$T>::try_from(-p).unwrap(), &m) }; // FIXME: implement ModInt for <T,T> and prevent this ugly call
                let q = if q >= 0 { <$T>::try_from(q).unwrap() % &m }
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
        }
    };
}

impl_lucasm_for_uprim!(u8);
impl_lucasm_for_uprim!(u16);
impl_lucasm_for_uprim!(u32);
impl_lucasm_for_uprim!(u64);
impl_lucasm_for_uprim!(u128);

impl LucasMod for BigUint {
    fn lucasm(p: isize, q: isize, m: Self, n: Self) -> (Self, Self) {
        let p = if p >= 0 { BigUint::from_isize(p).unwrap() % &m }
                        else { ModInt::<&BigUint, &BigUint>::negm(&BigUint::from_isize(-p).unwrap(), &m) };
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
}

impl<T: Integer + NumRef + Clone + FromPrimitive + LucasMod> PrimalityUtils for T
where for<'r> &'r T: RefNum<T> + std::ops::Shr<usize, Output = T> + ModInt<&'r T, &'r T, Output = T>
{
    fn is_prp(&self, base: Self) -> bool {
        let tm1 = self - Self::one();
        base.powm(&tm1, self).is_one()
    }

    fn is_sprp(&self, base: T) -> bool {
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

    fn is_lprp(&self, p: isize, q: isize) -> bool {
        // TODO: failed for p=1, self=25 or 125
        if self.is_even() { return false; }
        let d = p*p - 4*q;
        let d = if d > 0 { Self::from_isize(d).unwrap() }
                   else { ModInt::<&Self, &Self>::negm(&Self::from_isize(d).unwrap(), self) };
        let delta = match d.jacobi(self) {
            0 => self.clone(), -1 => self + Self::one(), 1 => self - Self::one(), _ => panic!("")
        };
        let (u, _) = LucasMod::lucasm(p, q, self.clone(), delta);
        u.is_zero()
    }

    fn is_slprp(&self, p: isize, q: isize) -> bool {
        // TODOL failed for p=1, self=25, 125 or 625
        if self.is_even() { return false; }
        let d = p*p - 4*q;
        let d = if d > 0 { Self::from_isize(d).unwrap() }
                   else { ModInt::<&Self, &Self>::negm(&Self::from_isize(d).unwrap(), self) };
        let delta = match d.jacobi(self) {
            0 => self.clone(), -1 => self + Self::one(), 1 => self - Self::one(), _ => panic!("")
        };

        let shift = delta.fac2();
        let base = &delta >> shift;

        let (ud, vd) = LucasMod::lucasm(p, q, self.clone(), base.clone());
        if ud.is_zero() || vd.is_zero() { return true; }

        // do first iteration to remove the sign on Q
        if shift == 0 { return false; }
        let qk = Self::from_isize(q.abs()).unwrap().powm(&base, self);
        let two = Self::from_u8(2).unwrap();
        // V(2k) = V(k)² - 2Q^k
        let mut vd = if q >= 0 { vd.mulm(&vd,self).subm(&two.mulm(&qk,self),self) }
                        else { vd.mulm(&vd,self).addm(&two.mulm(&qk,self),self) };
        let mut qk = qk.mulm(&qk,self);

        for _ in 1..shift {
            // V(2k) = V(k)² - 2Q^k
            vd = vd.mulm(&vd,self).subm(&two.mulm(&qk,self),self);
            if vd.is_zero() { return true }
            qk = qk.mulm(&qk,self);
        }
        false
    }

    fn is_eslprp(&self, p: isize) -> bool {
        unimplemented!()
    }

    fn is_slsprp(&self) -> bool {
        // TODO: also add perfect square check here
        unimplemented!()
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
            let (uk, _) = LucasMod::lucasm(3, -1, m as u64, n as u64);
            assert_eq!(uk, p3qm1[n] % (m as u64));
            let (uk, _) = LucasMod::lucasm(3, -1, BigUint::from(m), BigUint::from(n));
            assert_eq!(uk, BigUint::from(p3qm1[n] % (m as u64)));
        }
    }

    #[test]
    fn lucas_prp_test() {
        assert!(19u8.is_lprp(3, -1));

        // least lucas pseudo primes for Q=-1
        let plimit: [u16; 5] = [323, 35, 119, 9, 9];
        for (i, l) in plimit.iter().cloned().enumerate() {
            let p = (i + 1) as isize;
            assert!(l.is_lprp(p, -1));

            // test four random numbers under limit
            for _ in 0..10 {
                let n = random::<u16>() % l;
                if n < 2 || n.is_sprp(2) { continue } // skip real primes
                // assert!(!n.is_lprp(p, -1), "lucas prp test on {} with p = {}", n, p);
            }
        }

        // least strong lucas pseudoprimes for Q=-1
        let plimit: [u16; 3] = [4181, 169, 119];
        for (i, l) in plimit.iter().cloned().enumerate() {
            let p = (i + 1) as isize;
            assert!(l.is_slprp(p, -1));

            // test random numbers under limit
            for _ in 0..100 {
                let n = random::<u16>() % l;
                if n < 2 || (n.is_sprp(2) && n.is_sprp(2)) { continue } // skip real primes
                assert!(!n.is_slprp(p, -1), "strong lucas prp test on {} with p = {}", n, p);
            }
        }
    }
}
