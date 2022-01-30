//! Standalone number theoretic functions that can be used without prime cache

use crate::buffer::{NaiveBuffer, PrimeBufferExt};
use crate::factor::{pollard_rho, squfof, trial_division};
use crate::tables::{MOEBIUS_ODD, SMALL_PRIMES};
use crate::traits::{
    FactorizationConfig, Primality, PrimalityTestConfig, PrimalityUtils,
};
use crate::primality::{PrimalityBase, PrimalityRefBase};
use rand::random;
use std::collections::BTreeMap;
use std::convert::TryFrom;

#[cfg(feature = "big-table")]
use crate::tables::{MILLER_RABIN_BASE32, MILLER_RABIN_BASE64};

/// This function does fast primality test on a u64 integer. It's based on
/// deterministic Miller-rabin tests. if target is larger than 2^64 or more
/// controlled primality tests are desired, please use [is_prime()]
#[cfg(not(feature = "big-table"))]
pub fn is_prime64(target: u64) -> bool {
    // shortcuts
    if target < 1 {
        return false;
    }
    if target & 1 == 0 {
        return target == 2;
    }

    // first find in the prime list
    if let Ok(u) = u8::try_from(target) {
        return SMALL_PRIMES.binary_search(&u).is_ok();
    }

    // Then do a deterministic Miller-rabin test
    // The collection of witnesses are from http://miller-rabin.appspot.com/
    if let Ok(u) = u16::try_from(target) {
        // 2, 3 for u16 range
        return u.is_sprp(2) && u.is_sprp(3);
    }
    if let Ok(u) = u32::try_from(target) {
        // 2, 7, 61 for u32 range
        return u.is_sprp(2) && u.is_sprp(7) && u.is_sprp(61);
    }

    // 2, 325, 9375, 28178, 450775, 9780504, 1795265022 for u64 range
    const WITNESS64: [u64; 7] = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];
    WITNESS64.iter().all(|&x| target.is_sprp(x))
}

/// This function does very fast primality test on a u64 integer is a prime number. It's based on
/// deterministic Miller-rabin tests with hashing. if target is larger than 2^64 or more controlled
/// primality tests are desired, please use is_prime() or PrimeBuffer::is_prime()
#[cfg(feature = "big-table")]
pub fn is_prime64(target: u64) -> bool {
    if target < 8167 {
        return SMALL_PRIMES.binary_search(&(target as u16)).is_ok();
    }

    // 32bit test
    const MAGIC: u32 = 0xAD625B89;
    if let Ok(u) = u32::try_from(target) {
        let base = u.wrapping_mul(MAGIC) >> 24;
        return u.is_sprp(MILLER_RABIN_BASE32[base as usize]);
    }

    // 49bit test
    if !target.is_sprp(2) {
        return false;
    }
    let u = target as u32; // truncate
    let base = u.wrapping_mul(MAGIC) >> 18;
    if !target.is_sprp(MILLER_RABIN_BASE64[base as usize] as u64) {
        return false;
    }
    if target < (1u64 << 49) {
        return true;
    }

    // 64bit test
    const SECOND_BASES: [u64; 8] = [15, 135, 13, 60, 15, 117, 65, 29];
    let base = base >> 13;
    target.is_sprp(SECOND_BASES[base as usize])
}

pub fn factors64(target: u64) -> BTreeMap<u64, usize> {
    // TODO: improve factorization performance
    // REF: https://mathoverflow.net/questions/114018/fastest-way-to-factor-integers-260
    //      https://hal.inria.fr/inria-00188645v3/document
    //      https://github.com/coreutils/coreutils/blob/master/src/factor.c
    //      https://github.com/uutils/coreutils/blob/master/src/uu/factor/src/cli.rs
    //      https://github.com/elmomoilanen/prime-factorization
    //      https://github.com/radii/msieve
    let f2 = target.trailing_zeros(); // quick check on factors of 2
    if f2 == 0 && is_prime64(target) {
        let mut result = BTreeMap::new();
        result.insert(target, 1);
        return result;
    }

    // trial division using primes in the table
    let target = target >> f2;
    let piter = SMALL_PRIMES.iter().skip(1).map(|&p| p as u64); // skip 2
    let (mut result, factored) = trial_division(piter, target, None);
    if f2 > 0 {
        result.insert(2, f2 as usize);
    } // add back 2
    let residual = match factored {
        Ok(res) => {
            if res != 1 {
                result.insert(res, 1);
            }
            return result;
        }
        Err(res) => res,
    };

    // then try pollard's rho and SQUFOF methods util fully factored
    let mut todo = vec![residual];
    const SQUFOF_MULTIPLIERS: [u16; 16] = [
        1,
        3,
        5,
        7,
        11,
        3 * 5,
        3 * 7,
        3 * 11,
        5 * 7,
        5 * 11,
        7 * 11,
        3 * 5 * 7,
        3 * 5 * 11,
        3 * 7 * 11,
        5 * 7 * 11,
        3 * 5 * 7 * 11,
    ];
    while let Some(target) = todo.pop() {
        if is_prime64(target) {
            *result.entry(target).or_insert(0) += 1;
        } else {
            let mut i = 1usize;
            let divisor = loop {
                // try SQUFOF after 4 failed pollard rho trials
                if i % 5 == 0 && (i / 5) < SQUFOF_MULTIPLIERS.len() {
                    // TODO: check if the residual is a sqaure number before SQUFOF
                    if let Some(p) = squfof(&target, SQUFOF_MULTIPLIERS[i / 5] as u64) {
                        break p;
                    }
                } else {
                    let start = random::<u64>() % target;
                    let offset = random::<u64>() % target;
                    if let Some(p) = pollard_rho(&target, start, offset) {
                        break p;
                    }
                }
                i += 1;
            };
            todo.push(divisor);
            todo.push(target / divisor);
        }
    }
    result
}

/// This function re-exports [crate::buffer::PrimeBufferExt::is_prime()] with a default buffer distance
pub fn is_prime<T: PrimalityBase>(target: &T, config: Option<PrimalityTestConfig>) -> Primality
where
    for<'r> &'r T: PrimalityRefBase<T>,
{
    NaiveBuffer::new().is_prime(target, config)
}

/// This function re-exports [crate::buffer::PrimeBufferExt::factors()] with a default buffer instance
pub fn factors<T: PrimalityBase>(
    target: T,
    config: Option<FactorizationConfig>,
) -> Result<BTreeMap<T, usize>, Vec<T>>
where
    for<'r> &'r T: PrimalityRefBase<T>,
{
    NaiveBuffer::new().factors(target, config)
}

/// This function re-exports [NaiveBuffer::primorial()]
pub fn primorial<T: PrimalityBase + std::iter::Product>(n: usize) -> T {
    NaiveBuffer::new().primorial(n)
}

/// This function calculate the Möbius function of the input integer
/// It will panic if the factorization failed. If the input integer is
/// very hard to factorize, it's better to use the `factors` function to
/// control how the factorization is done.
pub fn moebius_mu<T: PrimalityBase>(target: &T) -> i8
where
    for<'r> &'r T: PrimalityRefBase<T>,
{
    // remove factor 2
    if target.is_even() {
        let two = T::one() + T::one();
        let four = &two + &two;
        if (target % four).is_zero() {
            return 0;
        } else {
            return -moebius_mu(&(target / &two));
        }
    }

    // look up tables when input is smaller than 256
    if let Some(v) = (target - T::one()).to_u8() {
        let m = MOEBIUS_ODD[(v >> 6) as usize];
        let m = m & (3 << (v & 63));
        let m = m >> (v & 63);
        return m as i8 - 1;
    }

    // short cut for common primes
    let three_sq = T::from_u8(9).unwrap();
    let five_sq = T::from_u8(25).unwrap();
    let seven_sq = T::from_u8(49).unwrap();
    if (target % three_sq).is_zero()
        || (target % five_sq).is_zero()
        || (target % seven_sq).is_zero()
    {
        return 0;
    }

    // then try complete factorization
    match factors(target.clone(), None) {
        Ok(result) => {
            for exp in result.values() {
                if exp > &1 {
                    return 0;
                }
            }
            return if result.len() % 2 == 0 { 1 } else { -1 };
        }
        Err(_) => {
            panic!("Failed to factor the integer!");
        }
    }
}

/// This function tests if the integer doesn't have any square number factor.
/// It will panic if the factorization failed.
pub fn is_square_free<T: PrimalityBase>(target: &T) -> bool
where
    for<'r> &'r T: PrimalityRefBase<T>,
{
    moebius_mu(target) != 0
}

// TODO: More functions
// REF: http://www.numbertheory.org/gnubc/bc_programs.html
// REF: https://github.com/TilmanNeumann/java-math-library
// - is_safe_prime: is safe prime with Sophie Germain definition
// - is_smooth: checks if the smoothness bound is at least b
// - euler_phi: Euler's totient function
// - jordan_tot: Jordan's totient function
// - moebius_mu: Möbius mu function
// Others include Louiville function, Mangoldt function, Dedekind psi function, Dickman rho function, etc..
//
// These function might be implemented in PrimeBuffer, ref http://flintlib.org/doc/ulong_extras.html#prime-number-generation-and-counting
// - prime_pi
// - prime_pi_bounds
// - nth_prime
// - nth_prime_bounds
// - next_prime

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{prelude::SliceRandom, random};
    use std::iter::FromIterator;

    #[test]
    fn is_prime64_test() {
        // test small primes
        for x in 2..100 {
            assert_eq!(SMALL_PRIMES.contains(&x), is_prime64(x as u64));
        }

        // some large primes
        assert!(is_prime64(6469693333));
        assert!(is_prime64(13756265695458089029));
        assert!(is_prime64(13496181268022124907));
        assert!(is_prime64(10953742525620032441));
        assert!(is_prime64(17908251027575790097));

        // primes from examples in Bradley Berg's hash method
        assert!(is_prime64(480194653));
        assert!(!is_prime64(20074069));
        assert!(is_prime64(8718775377449));
        assert!(is_prime64(3315293452192821991));
        assert!(!is_prime64(8651776913431));
        assert!(!is_prime64(1152965996591997761));

        // ensure no factor for 100 random primes
        let mut rng = rand::thread_rng();
        for _ in 0..100 {
            let x = random();
            if !is_prime64(x) {
                continue;
            }
            assert_ne!(x % (*SMALL_PRIMES.choose(&mut rng).unwrap() as u64), 0);
        }

        // create random composites
        for _ in 0..100 {
            let x = random::<u32>() as u64;
            let y = random::<u32>() as u64;
            assert!(!is_prime64(x * y));
        }
    }

    #[test]
    fn factors64_test() {
        // some known cases
        let fac123456789 = BTreeMap::from_iter([(3, 2), (3803, 1), (3607, 1)]);
        let fac = factors64(123456789);
        assert_eq!(fac, fac123456789);

        let fac1_17 = BTreeMap::from_iter([(2071723, 1), (5363222357, 1)]);
        let fac = factors64(11111111111111111);
        assert_eq!(fac, fac1_17);

        // 100 random factorization tests
        for _ in 0..100 {
            let x = random();
            let fac = factors64(x);
            let mut prod = 1;
            for (p, exp) in fac {
                assert!(
                    is_prime64(p),
                    "factorization result should have prime factors! (get {})",
                    p
                );
                prod *= p.pow(exp as u32);
            }
            assert_eq!(x, prod, "factorization check failed! ({} != {})", x, prod);
        }
    }

    #[test]
    fn moebius_mu_test() {
        // test small examples
        let mu20: [i8; 20] = [1,-1,-1,0,-1,1,-1,0,0,1,-1,0,-1,1,1,0,-1,0,-1,0];
        for i in 0..20 {
            assert_eq!(moebius_mu(&(i+1)), mu20[i], "moebius on {}", i);
        }

        // some square numbers
        assert_eq!(moebius_mu(&1024u32), 0);
        assert_eq!(moebius_mu(&(8081u32 * 8081)), 0);

        // sphenic numbers
        let sphenic3: [u8; 20] = [30, 42, 66, 70, 78, 102, 105, 110, 114, 130, 138, 154, 165, 170, 174, 182, 186, 190, 195, 222]; // OEIS A007304
        for i in 0..20 {
            assert_eq!(moebius_mu(&sphenic3[i]), -1i8, "moebius on {}", sphenic3[i]);
        }
        let sphenic5: [u16; 23] = [2310, 2730, 3570, 3990, 4290, 4830, 5610, 6006, 6090, 6270, 6510, 6630, 7410, 7590, 7770, 7854, 8610, 8778, 8970, 9030, 9282, 9570, 9690]; // OEIS A046387
        for i in 0..20 {
            assert_eq!(moebius_mu(&sphenic5[i]), -1i8, "moebius on {}", sphenic5[i]);
        }
    }
}
