use crate::factor::{pollard_rho, trial_division};
use crate::nt_funcs::{factors64, is_prime64};
use crate::primality::LucasUtils;
use crate::tables::{SMALL_PRIMES, SMALL_PRIMES_NEXT};
use crate::traits::{
    BitTest, FactorizationConfig, ModInt, Primality, PrimalityTestConfig, PrimalityUtils,
    PrimeBuffer,
};
use bitvec::bitvec;
use num_bigint::BigUint; // TODO (v0.1): make the dependency for this optional
use num_integer::{Integer, Roots};
use num_traits::{FromPrimitive, NumRef, RefNum, ToPrimitive};
use rand::{random, seq::IteratorRandom};
use std::{collections::BTreeMap, convert::TryInto};

pub trait PrimalityBase:
    Integer + Roots + NumRef + Clone + FromPrimitive + ToPrimitive + LucasUtils + BitTest
{
}
impl<T: Integer + Roots + NumRef + Clone + FromPrimitive + ToPrimitive + LucasUtils + BitTest>
    PrimalityBase for T
{
}
pub trait PrimalityRefBase<Base>:
    RefNum<Base>
    + std::ops::Shr<usize, Output = Base>
    + for<'r> ModInt<&'r Base, &'r Base, Output = Base>
{
}
impl<T, Base> PrimalityRefBase<Base> for T where
    T: RefNum<Base>
        + std::ops::Shr<usize, Output = Base>
        + for<'r> ModInt<&'r Base, &'r Base, Output = Base>
{
}

pub trait PrimeBufferExt: for<'a> PrimeBuffer<'a> {
    /// Test if an integer is a prime, the config will take effect only if the target is larger
    /// than 2^64.
    fn is_prime<T: PrimalityBase>(
        &self,
        target: &T,
        config: Option<PrimalityTestConfig>,
    ) -> Primality
    where
        for<'r> &'r T: PrimalityRefBase<T>,
    {
        // shortcuts
        if target.is_even() {
            return if target == &T::from_u8(2u8).unwrap() {
                Primality::Yes
            } else {
                Primality::No
            };
        }

        // do deterministic test if target is under 2^64
        if let Some(x) = target.to_u64() {
            return match is_prime64(x) {
                true => Primality::Yes,
                false => Primality::No,
            };
        }

        let config = config.unwrap_or(PrimalityTestConfig::default());
        let mut probability = 1.;

        // miller-rabin test
        let mut witness_list: Vec<u64> = Vec::new();
        if config.sprp_trials > 0 {
            witness_list.extend(self.iter().take(config.sprp_trials));
            probability *= 1. - 0.25f32.powi(config.sprp_trials.try_into().unwrap());
        }
        if config.sprp_random_trials > 0 {
            let mut rng = rand::thread_rng();
            witness_list.extend(
                self.iter()
                    .choose_multiple(&mut rng, config.sprp_random_trials),
            );
            probability *= 1. - 0.25f32.powi(config.sprp_random_trials.try_into().unwrap());
        }
        if !witness_list
            .iter()
            .all(|&x| target.is_sprp(T::from_u64(x).unwrap()))
        {
            return Primality::No;
        }

        // lucas probable prime test
        if config.slprp_test {
            probability *= 1. - 4f32 / 15f32;
            if !target.is_slprp(None, None) {
                return Primality::No;
            }
        }
        if config.eslprp_test {
            probability *= 1. - 4f32 / 15f32;
            if !target.is_eslprp(None) {
                return Primality::No;
            }
        }

        Primality::Probable(probability)
    }

    /// Return list of found factors if not fully factored
    fn factors<T: PrimalityBase>(
        &self,
        target: T,
        config: Option<FactorizationConfig>,
    ) -> Result<BTreeMap<T, usize>, Vec<T>>
    where
        for<'r> &'r T: PrimalityRefBase<T>,
    {
        // shortcut if the target is in u64 range
        if let Some(x) = target.to_u64() {
            return Ok(factors64(x)
                .iter()
                .map(|(&k, &v)| (T::from_u64(k).unwrap(), v))
                .collect());
        }
        let config = config.unwrap_or(FactorizationConfig::default());

        // test the existing primes
        let (result, factored) = trial_division(self.iter().cloned(), target, config.td_limit);
        let mut result: BTreeMap<T, usize> = result
            .into_iter()
            .map(|(k, v)| (T::from_u64(k).unwrap(), v))
            .collect();

        // find factors by dividing
        let mut successful = true;
        let mut config = config;
        config.td_limit = Some(0); // disable trial division when finding divisor
        match factored {
            Ok(res) => {
                if !res.is_one() {
                    result.insert(res, 1);
                }
            }
            Err(res) => {
                successful = true;
                let mut todo = vec![res];
                while let Some(target) = todo.pop() {
                    if self
                        .is_prime(&target, Some(config.prime_test_config))
                        .probably()
                    {
                        *result.entry(target).or_insert(0) += 1;
                    } else {
                        if let Some(divisor) = self.divisor(&target, &mut config) {
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

        if successful {
            Ok(result)
        } else {
            Err(result
                .into_iter()
                .flat_map(|(f, n)| std::iter::repeat(f).take(n))
                .collect())
        }
    }

    /// Return a proper divisor of target (randomly), even works for very large numbers
    /// Return None if no factor is found (this method will not do a primality check)
    fn divisor<T: PrimalityBase>(&self, target: &T, config: &mut FactorizationConfig) -> Option<T>
    where
        for<'r> &'r T: PrimalityRefBase<T>,
    {
        if matches!(config.td_limit, Some(0)) {
            // try to get a factor by trial division
            let tsqrt: T = Roots::sqrt(target) + T::one();
            let limit = if let Some(l) = config.td_limit {
                tsqrt.clone().min(T::from_u64(l).unwrap())
            } else {
                tsqrt.clone()
            };

            for p in self.iter().map(|p| T::from_u64(*p).unwrap()) {
                if &p > &tsqrt {
                    return None; // the number is a prime
                }
                if &p > &limit {
                    break;
                }
                if target.is_multiple_of(&p) {
                    return Some(p);
                }
            }
        }

        // try to get a factor using pollard_rho with 4x4 trials
        let below64 = target.to_u64().is_some();
        while config.rho_trials > 0 {
            let (start, offset) = if below64 {
                (
                    T::from_u8(random::<u8>()).unwrap() % target,
                    T::from_u8(random::<u8>()).unwrap() % target,
                )
            } else {
                (
                    T::from_u64(random::<u64>()).unwrap() % target,
                    T::from_u64(random::<u64>()).unwrap() % target,
                )
            };
            config.rho_trials -= 1;
            if let Some(p) = pollard_rho(target, start, offset) {
                return Some(p);
            }
        }

        None
    }
}

impl<T> PrimeBufferExt for T where for<'a> T: PrimeBuffer<'a> {}

/// NaiveBuffer implements a very simple Sieve of Eratosthenes
pub struct NaiveBuffer {
    list: Vec<u64>, // list of found prime numbers
    current: u64, // all primes smaller than this value has to be in the prime list, should be an odd number
}

impl NaiveBuffer {
    #[inline]
    pub fn new() -> Self {
        let list = SMALL_PRIMES.iter().map(|&p| p as u64).collect();
        NaiveBuffer {
            list,
            current: SMALL_PRIMES_NEXT,
        }
    }
}

impl<'a> PrimeBuffer<'a> for NaiveBuffer {
    type PrimeIter = std::slice::Iter<'a, u64>;

    fn contains(&self, num: u64) -> bool {
        self.list.binary_search(&num).is_ok()
    }

    fn clear(&mut self) {
        self.list.truncate(54); // reserve small primes below 255
        self.list.shrink_to_fit();
        self.current = 255;
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
        if odd_limit <= current {
            return;
        }

        // create sieve and filter with existing primes
        let mut sieve = bitvec![0; ((odd_limit - current) / 2) as usize];
        for p in self.list.iter().skip(1) {
            // skip pre-filtered 2
            let start = if p * p < current {
                p * ((current / p) | 1) // start from an odd factor
            } else {
                p * p
            };
            for multi in (start..odd_limit).step_by(2 * (*p as usize)) {
                if multi >= current {
                    sieve.set(((multi - current) / 2) as usize, true);
                }
            }
        }

        // sieve with new primes
        for p in (current..Roots::sqrt(&odd_limit) + 1).step_by(2) {
            for multi in (p * p..odd_limit).step_by(2 * (p as usize)) {
                if multi >= current {
                    sieve.set(((multi - current) / 2) as usize, true);
                }
            }
        }

        // sort the sieve
        self.list
            .extend(sieve.iter_zeros().map(|x| (x as u64) * 2 + current));
        self.current = odd_limit;
    }
}

impl NaiveBuffer {
    // FIXME: These two functions could be implemented in the trait, but only after
    // RFC 2071 and https://github.com/cramertj/impl-trait-goals/issues/3

    /// Returns all primes **below** limit. The primes are sorted.
    pub fn primes(&mut self, limit: u64) -> std::iter::Take<<Self as PrimeBuffer>::PrimeIter> {
        self.reserve(limit);
        let position = match self.list.binary_search(&limit) {
            Ok(p) => p,
            Err(p) => p,
        };
        return self.list.iter().take(position);
    }

    /// Returns primes of certain amount counting from 2. The primes are sorted.
    pub fn nprimes(&mut self, count: usize) -> std::iter::Take<<Self as PrimeBuffer>::PrimeIter> {
        loop {
            // TODO: use a more accurate function to estimate the upper/lower bound of prime number function pi(.)
            self.reserve(self.current * (count as u64) / (self.list.len() as u64));
            if self.list.len() >= count {
                break self.list.iter().take(count);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::random;
    use std::iter::FromIterator;

    const PRIME50: [u64; 15] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
    const PRIME100: [u64; 25] = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97,
    ];

    #[test]
    fn prime_generation_test() {
        let mut pb = NaiveBuffer::new();
        assert_eq!(pb.primes(50).cloned().collect::<Vec<_>>(), PRIME50);
        assert_eq!(pb.primes(100).cloned().collect::<Vec<_>>(), PRIME100);
    }

    #[test]
    fn pb_is_prime_test() {
        // test for is_prime
        let pb = NaiveBuffer::new();

        // some mersenne numbers
        assert!(matches!(
            pb.is_prime(&BigUint::from(2u32.pow(19) - 1), None),
            Primality::Yes
        ));
        assert!(matches!(
            pb.is_prime(&BigUint::from(2u32.pow(23) - 1), None),
            Primality::No
        ));
        let m89 = BigUint::from(2u8).pow(89) - 1u8;
        assert!(matches!(pb.is_prime(&m89, None), Primality::Probable(_)));

        // test against small prime assertion
        for _ in 0..100 {
            let target = random::<u64>();
            assert_eq!(
                !is_prime64(target),
                matches!(
                    pb.is_prime(&BigUint::from(target), Some(PrimalityTestConfig::bpsw())),
                    Primality::No
                )
            );
        }
    }

    #[test]
    fn pb_factors_test() {
        let pb = NaiveBuffer::new();

        let m131 = BigUint::from(2u8).pow(131) - 1u8; // m131/263 is a large prime
        let fac = pb.factors(m131, None);
        assert!(matches!(fac, Ok(f) if f.len() == 2));
    }
}
