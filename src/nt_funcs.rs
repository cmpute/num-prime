//! This module contains standalone number theoretic functions that can be used without prime cache

use crate::traits::Primality;
use std::collections::BTreeMap;

/// This function do BSW primality test on the target, if target is large enought
/// For more controlled primality test, please use NaiveBuffer::is_prime or NaiveBuffer::is_bprime
pub fn is_prime<T>(target: T) -> Primality {
    // targeting FLINT functions
    unimplemented!()
}

/// This function performs integer factorization on the target.
/// For more controlled factorization, please use NaiveBuffer::factors or NaiveBuffer::bfactors
pub fn factors<T>(target: T) -> Result<BTreeMap<T, usize>, Vec<T>>{
    // targeting GNU factors program and FLINT functions
    unimplemented!()
}

pub fn primes(limit: u64) { unimplemented!() }
pub fn nprimes(count: usize) { unimplemented!() }
pub fn primorial<T>(n: u64) -> T { unimplemented!() }

// TODO: More functions
// - is_prime64: fast prime check under 2^64, use either few miller-rabin or hash miller-rabin (create a feature option)
// - factor64: fast factorization under 2^64, use is_prime64 and maybe other tricks
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
