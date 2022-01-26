/// NaiveBuffer implements a list of primes

use std::{collections::BTreeMap, convert::{TryInto, TryFrom}};
use bitvec::bitvec;
use num_traits::{One, ToPrimitive, FromPrimitive, NumRef, RefNum};
use num_bigint::BigUint; // TODO: make the dependency for this optional
use num_integer::{Integer, Roots};
use rand::{random, seq::IteratorRandom};
use crate::traits::{PrimalityUtils, PrimeBuffer};
use crate::factor::pollard_rho;

pub enum Primality {
    Yes, No,
    /// carrying the probability of the number being a prime
    Probable(f32)
}

/// Find factors by trial division. The target is guaranteed fully factored
/// only if bound() * bound() > target. The parameter limit sets the max prime to be tried aside from bound()
/// Return factors, and if the target is fully factored, return Ok(residual), otherwise return Err(residual)
fn factors_trial<I: Iterator<Item = u64>, T: Integer + Clone + Roots + NumRef + FromPrimitive>
(primes: I, target: T, limit: Option<u64>) -> (BTreeMap<u64, usize>, Result<T, T>)
where for<'r> &'r T: RefNum<T> {
    let mut residual = target.clone();
    let mut result = BTreeMap::new();
    let target_sqrt = num_integer::sqrt(target) + T::one();
    let limit = if let Some(l) = limit { target_sqrt.clone().min(T::from_u64(l).unwrap()) } else { target_sqrt.clone() };
    let mut factored = false;
    for (p, pt) in primes.map(|p| (p, T::from_u64(p).unwrap())) {
        if &pt > &target_sqrt { factored = true; }
        if &pt > &limit { break; }

        while residual.is_multiple_of(&pt) {
            residual = residual / &pt;
            *result.entry(p).or_insert(0) += 1;
        }
        if residual == T::one() {
            break
        }
    }

    (result, if factored { Ok(residual) } else { Err(residual) })
}

#[derive(Clone, Copy)]
pub struct PrimalityTestConfig {
    pub sprp_trials: usize, // number of SPRP test, starting from base 2 
    pub sprp_random_trials: usize, // number of SPRP test with random base
    pub slprp_test: bool, // TODO(v0.0.4): Implement BPSW test
    pub eslprp_test: bool
}

impl PrimalityTestConfig {
    pub fn default() -> Self {
        Self { sprp_trials: 2, sprp_random_trials: 2, slprp_test: false, eslprp_test: false }
    }

    /// Create a configuration for Baillie-PSW test (base 2 SPRP test + SLPRP test)
    pub fn bpsw() -> Self {
        Self { sprp_trials: 1, sprp_random_trials: 0, slprp_test: true, eslprp_test: false }
    }

    /// Create a configuration for PSW test (base 2 SPRP + Fibonacci test)
    pub fn psw() { todo!() }
}

#[derive(Clone, Copy)]
pub struct FactorizationConfig {
    /// config for test if a 
    pub prime_test_config: PrimalityTestConfig,

    /// prime limit of trial division, you also need to reserve the buffer if all primes under the limit are to be tested.
    /// None means using all available primes
    pub tf_limit: Option<u64>,

    /// number of trials with Pollard's rho method
    pub rho_trials: usize,

    /// number of trials with Pollard's rho method (Brent variant)
    pub brent_trials: usize,

    /// number of trials with Pollard's p-1 method
    pub pm1_trials: usize,

    /// number of trials with William's p+1 method
    pub pp1_trials: usize,
}

impl FactorizationConfig {
    pub fn default() -> Self {
        Self {
            prime_test_config: PrimalityTestConfig::default(),
            tf_limit: Some(1 << 14), rho_trials: 4,
            brent_trials: 0, pm1_trials: 0, pp1_trials: 0
        }
    }
}

pub trait PrimeBufferExt : for<'a> PrimeBuffer<'a> {
    /// Return whether target is a prime. It uses Miller-rabin tests so even works
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
        let witnesses = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];
        witnesses.iter().all(|&x| target.is_sprp(x))
    }

    /// Test if a big integer is a prime, this function would carry out a probability test if
    /// the target is not smaller than 2^64
    fn is_bprime(&self, target: &BigUint, config: Option<PrimalityTestConfig>) -> Primality {
        // shortcuts
        if target.is_even() {
            return if target == &BigUint::from(2u8) { Primality::Yes } else { Primality::No };
        }

        // do deterministic test if target is under 2^64
        if let Some(x) = target.to_u64() {
            return match self.is_prime(x) {
                true => Primality::Yes, false => Primality::No
            };
        }

        let config = config.unwrap_or(PrimalityTestConfig::default());
        let mut probability = 1.;

        // miller-rabin test
        let mut witness_list: Vec<u64> = Vec::new();
        if config.sprp_trials > 0 {
            witness_list.extend(self.iter().take(config.sprp_trials));
            probability *= 1. - 0.25_f32.powi(config.sprp_trials.try_into().unwrap());
        }
        if config.sprp_random_trials > 0 {
            let mut rng = rand::thread_rng();
            witness_list.extend(self.iter().choose_multiple(&mut rng, config.sprp_random_trials));
            probability *= 1. - 0.25_f32.powi(config.sprp_random_trials.try_into().unwrap());
        }
        if !witness_list.iter().all(|&x| target.is_sprp(BigUint::from(x))) {
            return Primality::No;
        }

        Primality::Probable(probability)
    }

    fn factors(&self, target: u64) -> BTreeMap<u64, usize> {        
        // TODO: improve factorization performance
        // REF: https://github.com/coreutils/coreutils/blob/master/src/factor.c
        //      https://github.com/uutils/coreutils/blob/master/src/uu/factor/src/cli.rs
        //      https://github.com/elmomoilanen/prime-factorization
        //      https://github.com/radii/msieve
        if self.is_prime(target) {
            let mut result = BTreeMap::new();
            result.insert(target, 1);
            return result;
        }

        const FACTOR_THRESHOLD: u64 = 1 << 14; // TODO: this threshold should be optimized for only u64 integers
        let (mut result, factored) = factors_trial(self.iter().cloned(), target, Some(FACTOR_THRESHOLD));
        
        match factored {
            Ok(res) => { result.insert(res, 1); },
            Err(res) => {
                let mut todo = vec![res];
                while let Some(target) = todo.pop() {
                    if self.is_prime(target) {
                        *result.entry(target).or_insert(0) += 1;
                    } else {
                        let divisor = loop {
                            let start = random::<u64>() % target;
                            let offset = random::<u64>() % target;
                            if let Some(p) = pollard_rho(&target, start, offset) {
                                break p
                            }
                        };
                        todo.push(divisor);
                        todo.push(res / divisor);
                    }
                }
            }
        };
        result
    }

    /// Return list of found factors if not fully factored
    // TODO: accept general integer as input (thus potentially support other bigint such as crypto-bigint)
    // REF: https://pypi.org/project/primefac/
    fn bfactors(&mut self, target: BigUint, config: Option<FactorizationConfig>) -> Result<BTreeMap<BigUint, usize>, Vec<BigUint>> {
        // shortcut if the target is in u64 range
        if let Some(x) = target.to_u64() {
            return Ok(self.factors(x).iter().map(|(&k, &v)| (BigUint::from(k), v)).collect());
        }
        let config = config.unwrap_or(FactorizationConfig::default());

        // test the existing primes
        let (result, factored) = factors_trial(self.iter().cloned(), target, config.tf_limit);
        let mut result: BTreeMap<BigUint, usize> = result.into_iter().map(|(k, v)| (BigUint::from(k), v)).collect();

        // find factors by dividing
        let mut successful = true;
        let mut config = config;
        config.tf_limit = Some(0); // disable TF when finding divisor
        match factored {
            Ok(res) => { result.insert(res, 1); },
            Err(res) => {
                successful = true;
                let mut todo = vec![res];
                while let Some(target) = todo.pop() {
                    if matches! (self.is_bprime(&target, Some(config.prime_test_config)), Primality::Yes | Primality::Probable(_)) {
                        *result.entry(target).or_insert(0) += 1;
                    } else {
                        let offset = random::<u64>() % &target;
                        if let Some(divisor) = self.bdivisor(&target, &mut config) {
                            todo.push(divisor.clone());
                            todo.push(target / divisor);
                        } else {
                            *result.entry(target).or_insert(0) += 1;
                            successful = false;
                        }
                    }
                }
            }
        };

        if successful { Ok(result) }
        else { Err(result.into_iter().flat_map(|(f, n)| std::iter::repeat(f).take(n)).collect()) }
    }

    /// Return a proper divisor of target (randomly), even works for very large numbers
    /// Return None if no factor is found (this method will not do a primality check)
    fn bdivisor(&self, target: &BigUint, config: &mut FactorizationConfig) -> Option<BigUint> {
        // try to get a factor by trial division
        // TODO (v0.1): skip the sqrt if tf_limit^2 < target
        let target_sqrt: BigUint = num_integer::sqrt(target.clone()) + BigUint::one();
        let limit = if let Some(l) = config.tf_limit { target_sqrt.clone().min(BigUint::from_u64(l).unwrap()) } else { target_sqrt.clone() };

        for p in self.iter().map(|p| BigUint::from_u64(*p).unwrap()) {
            if &p > &target_sqrt { return None; } // the number is a prime
            if &p > &limit { break; }
            if target.is_multiple_of(&p) { return Some(p); }
        }

        // try to get a factor using pollard_rho with 4x4 trials
        while config.rho_trials > 0 {
            let start = random::<u64>() % target;
            let offset = random::<u64>() % target;
            config.rho_trials -= 1;
            if let Some(p) = pollard_rho(target, start, offset) {
                return Some(p);
            }
        }

        None
    }
}

impl<T> PrimeBufferExt for T where for <'a> T: PrimeBuffer<'a> {}

pub struct NaiveBuffer {
    list: Vec<u64>, // list of found prime numbers
    current: u64 // all primes smaller than this value has to be in the prime list, should be an odd number
}

impl NaiveBuffer {
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
    // FIXME: These two functions could be implemented in the trait, but only after
    // RFC 2071 and https://github.com/cramertj/impl-trait-goals/issues/3

    /// Returns all primes **below** limit. The primes are sorted.
    fn primes(&mut self, limit: u64) -> std::iter::Take<<Self as PrimeBuffer>::PrimeIter> {
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

    const PRIME50: [u64; 15] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
    const PRIME100: [u64; 25] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];

    #[test]
    fn prime_generation_test(){
        let mut pb = NaiveBuffer::new();
        assert_eq!(pb.primes(50).cloned().collect::<Vec<_>>(), PRIME50);
        assert_eq!(pb.primes(100).cloned().collect::<Vec<_>>(), PRIME100);
    }
    
    #[test]
    fn prime_assertion_test() {
        let mut pb = NaiveBuffer::new();
        assert!(pb.is_prime(6469693333));
        for x in 2..100 {
            assert_eq!(PRIME100.contains(&x), pb.is_prime(x));
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
        let m89 = BigUint::from(2u8).pow(89) - 1u8;
        assert!(matches!(pb.is_bprime(&m89, None), Primality::Probable(_)));
    }

    #[test]
    fn factorization_test() {
        let mut pb = NaiveBuffer::new();
        let fac123456789 = BTreeMap::from_iter([(3, 2), (3803, 1), (3607, 1)]);
        let fac = pb.factors(123456789);
        assert_eq!(fac, fac123456789);

        let m131 = BigUint::from(2u8).pow(131) - 1u8; // m131/263 is a large prime
        let fac = pb.bfactors(m131, None);
        assert!(matches!(fac, Ok(f) if f.len() == 2));
    }
}

