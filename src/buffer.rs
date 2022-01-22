/// NaiveBuffer implements a list of primes

use std::collections::BTreeMap;
use bitvec::prelude::bitvec;
use std::convert::TryFrom;
use num_traits::{ToPrimitive, One, Pow};
use num_bigint::BigUint; // TODO: make the dependency for this optional
use num_integer::Integer;
use rand::{random, seq::IteratorRandom};
use crate::traits::PrimeArithmetic;

pub enum Primality {
    Yes, No,
    /// carrying the probability of the number being a prime
    Probable(f32)
}

pub trait PrimeBuffer<'a> {
    type PrimeIter: Iterator<Item = &'a u64>;

    // directly return an iterator of existing primes
    fn iter(&'a self) -> Self::PrimeIter;

    // generate primes until the upper bound is equal or larger than limit
    fn reserve(&mut self, limit: u64);

    // get the upper bound of primes in the list
    fn bound(&self) -> u64;

    // test if the number is in the buffer
    fn contains(&self, num: u64) -> bool;

    // clear the prime buffer to save memory
    fn clear(&mut self);
}

pub trait PrimeBufferExt : for<'a> PrimeBuffer<'a> {
    /// Return whether target is a prime. It uses Miller test so even works
    /// for very large numbers and it's very fast
    fn is_prime(&self, target: u64) -> bool {
        assert!(target > 1);

        // shortcuts
        if target.is_even() {
            return target == 2;
        }

        // first find in the prime list
        if target < self.bound() {
            return self.contains(target);
        }

        // Then do a deterministic Miller-rabin test
        // TODO: implement is_sprp for u8,u16,u32
        // The collection of witnesses are from http://miller-rabin.appspot.com/
        if let Ok(u) = u16::try_from(target) {
            // 2, 3 for u16 range
            return target.is_sprp(2) && target.is_sprp(3);
        }
        if let Ok(u) = u32::try_from(target) {
            // 2, 7, 61 for u32 range
            return target.is_sprp(2) && target.is_sprp(7) && target.is_sprp(61);
        }
        
        // 2, 325, 9375, 28178, 450775, 9780504, 1795265022 for u64 range
        let witnesses = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];
        witnesses.iter().all(|&x| target.is_sprp(x))
    }

    /// Test if a big integer is a prime, this function would carry out a probability test
    /// If `trials` is positive, then witness numbers are selected randomly, otherwise selecting from start
    /// TODO: change trials to a config struct for more detailed control
    fn is_bprime(&self, target: &BigUint, trials: Option<i32>) -> Primality {
        // shortcuts
        if target.is_even() {
            return if target == &BigUint::from(2u8) { Primality::Yes } else { Primality::No };
        }

        if let Some(x) = target.to_u64() {
            return match self.is_prime(x) {
                true => Primality::Yes, false => Primality::No
            };
        }

        // miller-rabin test
        // TODO: improve based on https://gmplib.org/manual/Prime-Testing-Algorithm
        //       or use https://en.wikipedia.org/wiki/Baillie%E2%80%93PSW_primality_test
        // TODO: always run 2-base first, and then random base
        let trials = trials.unwrap_or(4);
        let witness_list = if trials > 0 { 
            let mut rng = rand::thread_rng();
            self.iter().choose_multiple(&mut rng, trials as usize)
        } else {
            self.iter().take((-trials) as usize).collect()
        };
        match witness_list.iter().all(|&x| target.is_sprp(BigUint::from(*x))) {
            true => Primality::Probable(1. - 0.25_f32.powi(trials.abs())), false => Primality::No
        }
    }

    fn factors(&mut self, target: u64) -> BTreeMap<u64, usize> {        
        // TODO: improve factorization performance
        // REF: https://github.com/coreutils/coreutils/blob/master/src/factor.c
        //      https://pypi.org/project/primefac/
        //      https://github.com/uutils/coreutils/blob/master/src/uu/factor/src/cli.rs
        if self.is_prime(target) {
            let mut result = BTreeMap::new();
            result.insert(target, 1);
            return result;
        }

        const FACTOR_THRESHOLD: u64 = 1 << 28;
        if target < FACTOR_THRESHOLD {
            self.factors_naive(target)
        } else {
            self.factors_divide(target)
        }
    }

    /// Return list of found factors if not fully factored
    /// `trials` determines the maximum Pollard rho trials for each component
    /// TODO: accept general integer as input (thus potentially support other bigint such as crypto-bigint)
    fn bfactors(&mut self, target: &BigUint, trials: Option<i32>) -> Result<BTreeMap<BigUint, usize>, Vec<BigUint>> {
        // if the target is in u64 range
        if let Some(x) = target.to_u64() {
            return Ok(self.factors(x).iter().map(|(&k, &v)| (BigUint::from(k), v)).collect());
        }

        // test the existing primes
        let mut residual = target.clone();
        let mut trivial = BTreeMap::new();
        for &p in self.iter() {
            while residual.is_multiple_of(&BigUint::from(p)) {
                residual /= p;
                *trivial.entry(p).or_insert(0) += 1;
            }
            if residual == BigUint::one() {
                return Ok(trivial.iter().map(|(&k, &v)| (BigUint::from(k), v)).collect());
            }
        }

        // find factors by dividing
        let divided = self.bfactors_divide(&residual, trials.unwrap_or(4));

        // check if the number is fully factored
        let mut verify = BigUint::one();
        for (factor, exponent) in &trivial {
            verify *= factor.pow(exponent);
        }
        for (factor, exponent) in &divided {
            verify *= Pow::pow(factor, exponent);
        }

        // return results
        if &verify == target {
            let mut result = divided;
            for (factor, exponent) in trivial {
                *result.entry(BigUint::from(factor)).or_insert(0) += &exponent;
            }
            Ok(result)
        } else {
            Err(trivial.into_iter().flat_map(|(f, n)| std::iter::repeat(BigUint::from(f)).take(n)).collect())
        }
    }

    fn factors_naive(&mut self, target: u64) -> BTreeMap<u64, usize> {
        debug_assert!(!self.is_prime(target));

        let mut residual = target;
        let mut result = BTreeMap::new();
        let limit = num_integer::sqrt(target) + 1;
        self.reserve(limit);
        for &p in self.iter().take_while(|&p| p < &limit) {
            while residual % p == 0 {
                residual /= p;
                *result.entry(p).or_insert(0) += 1;
            }
            if residual == 1 {
                break
            }
        }

        if residual != 1 {
            result.insert(residual, 1);
        }
        result
    }

    /// Find the factors by dividing the target by a proper divider recursively
    fn factors_divide(&mut self, target: u64) -> BTreeMap<u64, usize> {
        debug_assert!(!self.is_prime(target));

        let d = self.divisor_rho(target);
        let mut f1 = self.factors(d);
        let f2 = self.factors(target / d);
        for (factor, exponent) in f2 {
            *f1.entry(factor).or_insert(0) += &exponent;
        }
        f1
    }

    /// Find the factors by dividing the target by a proper divider, if the target is prime
    /// or no divider is found, then an empty map is returned.
    ///
    /// Note: 
    /// We don't factorize probable prime since it will takes a long time.
    /// To factorize a probable prime, use bdivisor
    fn bfactors_divide(&self, target: &BigUint, trials: i32) -> BTreeMap<BigUint, usize> {
        if matches! (self.is_bprime(target, None), Primality::Yes | Primality::Probable(_)) {
            return BTreeMap::new();
        }

        match self.bdivisor_rho(target, trials) {
            Some(d) => {
                let mut f1 = self.bfactors_divide(&d, trials);
                if f1.len() == 0 { f1.insert(d.clone(), 1); } // add divisor if it's a prime
                let f2 = self.bfactors_divide(&(target / d), trials);
                for (factor, exponent) in f2 {
                    *f1.entry(factor).or_insert(0) += exponent;
                }
                f1
            },
            None => BTreeMap::new()
        }
    }

    /// Return a proper divisor of target (randomly), even works for very large numbers
    /// Return None if it's a prime
    fn divisor(&mut self, target: u64) -> Option<u64> {
        if self.is_prime(target) { return None }

        const DIVISOR_THRESHOLD: u64 = 1 << 36;
        if target < DIVISOR_THRESHOLD {
            Some(self.divisor_naive(target))
        } else {
            Some(self.divisor_rho(target))
        }
    }

    /// Return a proper divisor of target (randomly), even works for very large numbers
    /// Return None if it's a prime or no factor is found
    /// `trials` determine max Pollard rho trials
    fn bdivisor(&mut self, target: &BigUint, trials: Option<i32>) -> Option<BigUint> {
        // if the target is in u64 range
        if let Some(x) = target.to_u64() {
            return self.divisor(x).map(BigUint::from);
        }

        let primality = self.is_bprime(target, None);
        if let Primality::Yes = primality {
            return None;
        }

        // try to get a factor using pollard_rho with 4x4 trials
        self.bdivisor_rho(target, trials.unwrap_or(4))
    }

    // Get a factor by naive trials
    fn divisor_naive(&mut self, target: u64) -> u64 {
        debug_assert!(!self.is_prime(target));

        let limit = num_integer::sqrt(target) + 1;
        self.reserve(limit);
        *self.iter().take_while(|&p| p < &limit)
            .filter(|&x| target % *x == 0)
            .next().unwrap()
    }

    // Get a factor using pollard_rho
    fn divisor_rho(&self, target: u64) -> u64 {
        debug_assert!(!self.is_prime(target));
        loop {
            let offset = random::<u64>() % target;
            if let Some(p) = target.pollard_rho(offset, 4) {
                break p
            }
        }
    }
    fn bdivisor_rho(&self, target: &BigUint, trials: i32) -> Option<BigUint> {
        for _ in 0..trials {
            let offset = random::<u64>();
            if let Some(p) = target.pollard_rho(BigUint::from(offset), 4) {
                return Some(p)
            }
        }
        None
    }

    // Calculates the primorial number
    // fn primorial(&mut self, n: usize) -> BigUint {
    //     self.nprimes(n).cloned().map(BigUint::from).product()
    // }
}

impl<T> PrimeBufferExt for T where for <'a> T: PrimeBuffer<'a> {}

pub struct NaiveBuffer {
    list: Vec<u64>, // list of found prime numbers
    current: u64 // all primes smaller than this value has to be in the prime list, should be an odd number
}

// TODO: create static functions that can do primality test and factorization, might need a trait
// TODO: add a config struct for prime test and factorization. add a function to automatically select params

impl NaiveBuffer {
    // TODO: support indexing and iterating
    // TODO: make this struct a trait, and use our implementation as NaiveBuffer, and then optionally support primal-sieve?
    #[inline]
    pub fn new() -> Self {
        // store at least enough primes for miller test
        let list = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
        NaiveBuffer { list , current: 41 }
    }
}

impl<'a> PrimeBuffer<'a> for NaiveBuffer {
    type PrimeIter = std::slice::Iter<'a, u64>;

    fn contains(&self, num: u64) -> bool {
        self.list.binary_search(&num).is_ok()
    }

    fn clear(&mut self) {
        self.list.truncate(12); // reserve 2 ~ 37 for miller test
        self.list.shrink_to_fit();
        self.current = 41;
    }

    fn iter(&'a self) -> Self::PrimeIter {
        self.list.iter()
    }

    fn bound(&self) -> u64 {
        *self.list.last().unwrap()
    }

    fn reserve(&mut self, limit: u64) {
        let odd_limit = limit | 1; // make sure limit is odd
        let current = self.current; // prevent borrowing self
        debug_assert!(current % 2 == 1);

        // create sieve and filter with existing primes
        let mut sieve = bitvec![0; ((odd_limit - current) / 2) as usize];
        for p in self.list.iter().skip(1) { // skip pre-filtered 2
            let start = if p * p < current {
                p * ((current / p) | 1) // start from an odd factor
            } else {
                p * p
            };
            for multi in (start .. odd_limit).step_by(2 * (*p as usize)) {
                if multi >= current {
                    sieve.set(((multi - current) / 2) as usize, true);
                }
            }
        }

        // sieve with new primes
        for p in (current..num_integer::sqrt(odd_limit) + 1).step_by(2) {
            for multi in (p*p .. odd_limit).step_by(2 * (p as usize)) {
                if multi >= current {
                    sieve.set(((multi - current) / 2) as usize, true);
                }
            }
        }

        // sort the sieve
        self.list.extend(sieve.iter_zeros().map(|x| (x as u64) * 2 + current));
        self.current = odd_limit;
    }
}

impl NaiveBuffer {
    /// Returns all primes **below** limit. The primes are sorted.
    fn primes(&mut self, limit: u64) -> std::iter::Take<<Self as PrimeBuffer>::PrimeIter> {
        // TODO: how to return a take_while iterator? in that way we can eliminate the binary search
        // After solving this problem we can put these two methods into PrimeBufferExt trait
        self.reserve(limit);
        let position = match self.list.binary_search(&limit) {
            Ok(p) => p, Err(p) => p
        };
        return self.list.iter().take(position)
    }

    /// Returns primes of certain amount counting from 2. The primes are sorted.
    fn nprimes(&mut self, count: usize) -> std::iter::Take<<Self as PrimeBuffer>::PrimeIter> {
        loop {
            // TODO: use a more accurate function to estimate the upper/lower bound of prime number function pi(.)
            self.primes(self.current * (count as u64) / (self.list.len() as u64));
            if self.list.len() >= count {
                break self.list.iter().take(count)
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::iter::FromIterator;

    #[test]
    fn prime_generation_test(){
        let mut pb = NaiveBuffer::new();
        let prime50 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
        let prime100 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];
        assert_eq!(pb.primes(50).cloned().collect::<Vec<_>>(), prime50);
        assert_eq!(pb.primes(100).cloned().collect::<Vec<_>>(), prime100);
    }
    
    #[test]
    fn prime_assertion_test() {
        let mut pb = NaiveBuffer::new();
        assert!(pb.is_prime(6469693333));
        let prime100 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];
        for x in 2..100 {
            assert_eq!(prime100.contains(&x), pb.is_prime(x));
        }

        // test primes under 10000
        let gprimes = pb.primes(10000).cloned().collect::<Vec<_>>();
        for x in gprimes {
            assert!(pb.is_prime(x));
        }

        // test primes under 20000
        let gprimes = pb.primes(20000).cloned().collect::<Vec<_>>();
        for x in gprimes {
            assert!(pb.is_prime(x));
        }

        assert!(matches!(pb.is_bprime(&BigUint::from(2u32.pow(19) - 1), None), Primality::Yes));
        assert!(matches!(pb.is_bprime(&BigUint::from(2u32.pow(23) - 1), None), Primality::No));
        let m89 = BigUint::from(2u8).pow(89usize) - 1u8;
        assert!(matches!(pb.is_bprime(&m89, None), Primality::Probable(_)));
    }

    #[test]
    fn factorization_test() {
        let mut pb = NaiveBuffer::new();
        let fac123456789 = BTreeMap::from_iter([(3, 2), (3803, 1), (3607, 1)]);
        let fac = pb.factors(123456789);
        assert_eq!(fac, fac123456789);

        let m131 = BigUint::from(2u8).pow(131usize) - 1u8; // m131/263 is a large prime
        let fac = pb.bfactors(&m131, None);
        assert!(matches!(fac, Err(f) if f.len() > 0));
    }
}

