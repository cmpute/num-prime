/// PrimeBuffer implements a list of primes

use std::collections::BTreeMap;
use bitvec::prelude::bitvec;
use num_traits::{ToPrimitive, One, Pow};
use num_bigint::BigUint; // TODO: make the dependency for this optional
use num_integer::Integer;
use rand::{random, seq::IteratorRandom};
use crate::traits::PrimeArithmetic;

pub struct PrimeBuffer {
    list: Vec<u64>, // list of found prime numbers
    current: u64 // all primes smaller than this value has to be in the prime list, should be an odd number
}

pub enum Primality {
    Yes, No,
    /// carrying the probability of the number being a prime
    Probable(f32)
}

// TODO: create static functions that can do primality test and factorization, might need a trait
// TODO: add a config struct for prime test and factorization. add a function to automatically select params

impl PrimeBuffer { // TODO: support indexing and iterating
    // TODO: is the PrimeSet in rust primal / C++ primesieve more memory efficient?
    #[inline]
    pub fn new() -> Self {
        // store at least enough primes for miller test
        let list = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
        PrimeBuffer { list , current: 41 }
    }

    /// Return whether target is a prime. It uses Miller test so even works
    /// for very large numbers and it's very fast
    pub fn is_prime(&self, target: u64) -> bool {
        assert!(target > 1);

        // shortcuts
        if target.is_even() {
            return target == 2;
        }

        // first find in the prime list
        if target < self.current {
            return self.list.binary_search(&target).is_ok();
        }

        // Then do a deterministic Miller test
        // TODO: use a small set for test
        // REF: https://github.com/uutils/coreutils/blob/master/src/uu/factor/src/miller_rabin.rs#L28
        //      http://miller-rabin.appspot.com/
        // 2, 3 for u16
        // 2, 7, 61 for u32
        // 2, 325, 9375, 28178, 450775, 9780504, 1795265022 for u64
        let max_a = match target {  // https://oeis.org/A014233
            0..=2047 => 2,
            2048..=1373653 => 3,
            1373654..=25326001 => 5,
            25326002..=3251031751 => 7,
            3251031752..=2152302898747 => 11,
            2152302898748..=3474749660383 => 13,
            3474749660384..=341550071728321 => 17,
            341550071728322..=3825123056546413051 => 23,
            3825123056546413052.. => 37
        };

        self.list.iter()
            .take_while(|&x| *x <= max_a)
            .all(|&x| target.is_sprp(x))
    }

    /// Test if a big integer is a prime, this function would carry out a probability test
    /// If `trials` is positive, then witness numbers are selected randomly, otherwise selecting from start
    /// TODO: change trials to a config struct for more detailed control
    pub fn is_bprime(&self, target: &BigUint, trials: Option<i32>) -> Primality {
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
            self.list.iter().choose_multiple(&mut rng, trials as usize)
        } else {
            self.list.iter().take((-trials) as usize).collect()
        };
        match witness_list.iter().all(|&x| target.is_sprp(BigUint::from(*x))) {
            true => Primality::Probable(1. - 0.25_f32.powi(trials.abs())), false => Primality::No
        }
    }

    pub fn factors(&mut self, target: u64) -> BTreeMap<u64, usize> {        
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
    pub fn bfactors(&mut self, target: &BigUint, trials: Option<i32>) -> Result<BTreeMap<BigUint, usize>, Vec<BigUint>> {
        // if the target is in u64 range
        if let Some(x) = target.to_u64() {
            return Ok(self.factors(x).iter().map(|(&k, &v)| (BigUint::from(k), v)).collect());
        }

        // test the existing primes
        let mut residual = target.clone();
        let mut trivial = BTreeMap::new();
        for &p in &self.list {
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

    pub fn factors_naive(&mut self, target: u64) -> BTreeMap<u64, usize> {
        debug_assert!(!self.is_prime(target));

        let mut residual = target;
        let mut result = BTreeMap::new();
        for &p in self.primes(num_integer::sqrt(target) + 1) {
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
    pub fn factors_divide(&mut self, target: u64) -> BTreeMap<u64, usize> {
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
    fn bfactors_divide(&mut self, target: &BigUint, trials: i32) -> BTreeMap<BigUint, usize> {
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
    pub fn divisor(&mut self, target: u64) -> Option<u64> {
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
    pub fn bdivisor(&mut self, target: &BigUint, trials: Option<i32>) -> Option<BigUint> {
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
    pub fn divisor_naive(&mut self, target: u64) -> u64 {
        debug_assert!(!self.is_prime(target));
        *self.primes(num_integer::sqrt(target) + 1).iter()
            .filter(|&x| target % x == 0)
            .next().unwrap()
    }

    // Get a factor using pollard_rho
    pub fn divisor_rho(&self, target: u64) -> u64 {
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

    /// Returns all primes **below** limit. The primes are sorted.
    /// 
    /// # Reference
    /// https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n/3035188#3035188
    /// 
    /// # Note
    /// List primes works very slow for limit larger than 2^25 and won't work if limit > 2^30
    pub fn primes(&mut self, limit: u64) -> &[u64] {
        let odd_limit = limit | 1; // make sure limit is odd
        let current = self.current; // prevent borrowing self
        debug_assert!(current % 2 == 1);

        if odd_limit <= current {
            let position = match self.list.binary_search(&odd_limit) {
                Ok(p) => p, Err(p) => p
            };
            return &self.list[..position]
        }

        // create sieve and filter with existing primes
        // TODO: use wheel30 / wheel210 to speed up
        let mut sieve = bitvec![0; ((odd_limit - current) / 2) as usize];
        for p in self.list.iter().skip(1) { // skip pre-filtered 2
            let start = if p * p < current {
                p * ((current / p) | 1)
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
            if sieve[(p - current) as usize] {
                continue;
            }
            for multi in (p*p .. odd_limit).step_by(2 * (p as usize)) {
                if multi >= current {
                    sieve.set(((multi - current) / 2) as usize, true);
                }
            }
        }

        // sort the sieve
        self.list.extend(sieve.iter_zeros().map(|x| (x as u64) * 2 + current));
        self.current = odd_limit;

        if *self.list.last().unwrap() > limit {
            &self.list[..self.list.len()-1] // in this case odd_limit = limit + 1 is prime
        } else {
            &self.list[..]
        }
    }

    /// Returns primes of certain amount counting from 2. The primes are sorted.
    pub fn nprimes(&mut self, count: usize) -> &[u64] {
        loop {
            // TODO: use a more accurate function to estimate the upper/lower bound of prime number function pi(.)
            self.primes(self.current * (count as u64) / (self.list.len() as u64));
            if self.list.len() >= count {
                break &self.list[..count]
            }
        }
    }

    // Calculates the primorial number
    pub fn primorial(&mut self, n: usize) -> BigUint {
        self.nprimes(n).iter().map(|&p| BigUint::from(p)).product()
    }

    pub fn clear(&mut self) {
        self.list.truncate(12); // reserve 2 ~ 37 for miller test
        self.list.shrink_to_fit();
        self.current = 41;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::iter::FromIterator;

    #[test]
    fn prime_generation_test(){
        let mut pb = PrimeBuffer::new();
        let prime50 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
        let prime100 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];
        assert_eq!(pb.primes(50), prime50);
        assert_eq!(pb.primes(100), prime100);
    }
    
    #[test]
    fn prime_assertion_test() {
        let mut pb = PrimeBuffer::new();
        assert!(pb.is_prime(6469693333));
        let prime100 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];
        for x in 2..100 {
            assert_eq!(prime100.contains(&x), pb.is_prime(x));
        }

        let gprimes = pb.primes(10000).to_vec();
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
        let mut pb = PrimeBuffer::new();
        let fac123456789 = BTreeMap::from_iter([(3, 2), (3803, 1), (3607, 1)]);
        let fac = pb.factors(123456789);
        assert_eq!(fac, fac123456789);

        let m131 = BigUint::from(2u8).pow(131usize) - 1u8; // m131/263 is a large prime
        let fac = pb.bfactors(&m131, None);
        assert!(matches!(fac, Err(f) if f.len() > 0));
    }
}

