//! This module contains standalone number theoretic functions that can be used without prime cache
use crate::factor::pollard_rho;
use crate::tables::SMALL_PRIMES;
use crate::traits::{Primality, PrimalityUtils};
use rand::random;
use std::collections::BTreeMap;
use std::convert::TryFrom;

#[cfg(feature = "big-table")]
use crate::tables::{MILLER_RABIN_BASE32, MILLER_RABIN_BASE64};

/// This function does fast primality test on a u64 integer is a prime number. It's based on
/// deterministic Miller-rabin tests. if target is larger than 2^64 or more controlled primality
/// tests are desired, please use is_prime() or PrimeBuffer::is_prime()
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
    const WITNESSES: [u64; 7] = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];
    WITNESSES.iter().all(|&x| target.is_sprp(x))
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
    let mut result = BTreeMap::new();

    // TODO: improve factorization performance
    // REF: https://github.com/coreutils/coreutils/blob/master/src/factor.c
    //      https://github.com/uutils/coreutils/blob/master/src/uu/factor/src/cli.rs
    //      https://github.com/elmomoilanen/prime-factorization
    //      https://github.com/radii/msieve
    if is_prime64(target) {
        result.insert(target, 1);
        return result;
    }

    // trial division using primes in the table
    let mut residual = target;
    let mut result = BTreeMap::new();
    let target_sqrt = num_integer::sqrt(target) + 1;
    let mut factored = false;
    for p in SMALL_PRIMES {
        // TODO: remove factor 2 by counting trailing_zeros() and skip 2
        let p = p as u64;
        if p >= target_sqrt {
            factored = true;
            break;
        }

        while residual % p == 0 {
            residual = residual / p;
            *result.entry(p).or_insert(0) += 1;
        }
        if residual == 1 {
            factored = true;
            break;
        }
    }

    if factored {
        if residual > 1 {
            result.insert(residual, 1);
        }
        return result;
    }

    // then try pollard's rho method util fully factored
    let mut todo = vec![residual];
    while let Some(target) = todo.pop() {
        if is_prime64(target) {
            *result.entry(target).or_insert(0) += 1;
        } else {
            let divisor = loop {
                let start = random::<u64>() % target;
                let offset = random::<u64>() % target;
                if let Some(p) = pollard_rho(&target, start, offset) {
                    break p;
                }
            };
            todo.push(divisor);
            todo.push(target / divisor);
        }
    }
    result
}

/// This function do BSW primality test on the target
pub fn is_prime<T>(target: T) -> Primality {
    unimplemented!()
}

/// This function performs integer factorization on the target.
pub fn factors<T>(target: T) -> Result<BTreeMap<T, usize>, Vec<T>> {
    unimplemented!()
}

pub fn primes(limit: u64) {
    unimplemented!()
}
pub fn nprimes(count: usize) {
    unimplemented!()
}
pub fn primorial<T>(n: u64) -> T {
    unimplemented!()
}

// TODO: More functions
// - is_safe_prime: is safe prime with Sophie Germain definition
// - is_smooth: checks if the smoothness bound is at least b
// - euler_phi: Euler's totient function
// - jordan_tot: Jordan's totient function
// - moebius_mu: Möbius mu function
// Others include Louiville function, Mangoldt function, Dedekind psi function, etc..
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
    use rand::random;
    use std::iter::FromIterator;

    #[test]
    fn is_prime64_test() {
        // test for is_prime
        assert!(is_prime64(6469693333));

        // primes from examples in Bradley Berg's hash method
        assert!(is_prime64(480194653));
        assert!(!is_prime64(20074069));
        assert!(is_prime64(8718775377449));
        assert!(is_prime64(3315293452192821991));
        assert!(!is_prime64(8651776913431));
        assert!(!is_prime64(1152965996591997761));

        for x in 2..100 {
            assert_eq!(SMALL_PRIMES.contains(&x), is_prime64(x as u64));
        }
    }

    #[test]
    fn factors64_test() {
        let fac123456789 = BTreeMap::from_iter([(3, 2), (3803, 1), (3607, 1)]);
        let fac = factors64(123456789);
        assert_eq!(fac, fac123456789);
    }
}
