mod buffer;
mod traits;
mod integer;

pub use buffer::{PrimeBuffer, PrimeBufferExt, NaiveBuffer, Primality, PrimalityTestConfig, FactorizationConfig};
pub use traits::{NumberTheoretic, ModInt};
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
