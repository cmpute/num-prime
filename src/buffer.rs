//! Implementations and extensions for [PrimeBuffer], which represents a container of primes
//!
//! In `num-prime`, there is no global instance to store primes, the user has to generate
//! and store the primes themselves. The trait [PrimeBuffer] defines a unified interface
//! for a prime number container. Some methods that can take advantage of pre-generated
//! primes will be implemented in the [PrimeBufferExt] trait.
//!
//! We also provide [NaiveBuffer] as a simple implementation of [PrimeBuffer] without any
//! external dependencies. The performance of the [NaiveBuffer] will not be extremely optimized,
//! but it will be efficient enough for most applications.
//!

use crate::factor::{pollard_rho, trial_division};
use crate::nt_funcs::{factors64, is_prime64, nth_prime_bounds};
use crate::primality::{PrimalityBase, PrimalityRefBase};
use crate::tables::{SMALL_PRIMES, SMALL_PRIMES_NEXT};
use crate::traits::{
    FactorizationConfig, Primality, PrimalityTestConfig, PrimalityUtils, PrimeBuffer,
};
use bitvec::bitvec;
use num_integer::Roots;
use rand::{random, seq::IteratorRandom};
use std::{collections::BTreeMap, convert::TryInto};

/// Extension functions that can utilize pre-generated primes
pub trait PrimeBufferExt: for<'a> PrimeBuffer<'a> {
    /// Test if an integer is a prime. The `config` will take effect only if the target is larger
    /// than 2^64, otherwise [is_prime64] will be used.
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

    /// Factorize an integer. The `config` will take effect only if the target is larger
    /// than 2^64, otherwise [factors64] will be used.
    ///
    /// The factorization result will be returned as a map from primes to exponents. If the
    /// factorization failed, then a list of found factors will be returned.
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

    /// Return a proper divisor of target (randomly), even works for very large numbers.
    /// Return `None` if no factor is found.
    ///
    /// Note: this method will not do a primality check
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
    // FIXME: These functions could be implemented in the trait, but only after
    // RFC 2071 and https://github.com/cramertj/impl-trait-goals/issues/3

    /// This function calculates primorial function on n
    pub fn primorial<T: PrimalityBase + std::iter::Product>(&mut self, n: usize) -> T {
        self.nprimes(n).map(|&p| T::from_u64(p).unwrap()).product()
    }

    /// Returns all primes ≤ `limit`. The primes are sorted.
    pub fn primes(&mut self, limit: u64) -> std::iter::Take<<Self as PrimeBuffer>::PrimeIter> {
        self.reserve(limit);
        let position = match self.list.binary_search(&limit) {
            Ok(p) => p + 1,
            Err(p) => p,
        }; // into_ok_or_err()
        return self.list.iter().take(position);
    }

    /// Returns primes of certain amount counting from 2. The primes are sorted.
    pub fn nprimes(&mut self, count: usize) -> std::iter::Take<<Self as PrimeBuffer>::PrimeIter> {
        let (_, bound) = nth_prime_bounds(&(count as u64))
            .expect("Estimated size of the largest prime will be larger than u64 limit");
        self.reserve(bound);
        debug_assert!(self.list.len() >= count);
        self.list.iter().take(count)
    }

    /// Calculate and return the nth prime
    pub fn nth_prime(&mut self, n: usize) -> u64 {
        *self.nprimes(n).last().unwrap()
    }

    /// Calculate and return the prime pi function, i.e. number of primes ≤ `limit`.
    pub fn prime_pi(&mut self, limit: u64) -> usize {
        self.reserve(limit);

        // Count from the end of the list. Binary search is not used here since usually
        // this function is called without reserve() in advance.
        let fcount = self.list.iter().rev().take_while(|&p| p > &limit).count();
        self.list.len() - fcount
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "num-bigint")]
    use core::str::FromStr;
    #[cfg(feature = "num-bigint")]
    use num_bigint::BigUint;
    use rand::random;

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

        // test when limit itself is a prime
        pb.clear();
        assert_eq!(pb.primes(97).cloned().collect::<Vec<_>>(), PRIME100);

        pb.clear();
        assert_eq!(pb.nth_prime(10000), 104729);
        assert_eq!(pb.nth_prime(20000), 224737);
        assert_eq!(pb.nth_prime(10000), 104729); // use existing primes

        pb.clear();
        assert_eq!(pb.prime_pi(10000), 1229);
        assert_eq!(pb.prime_pi(20000), 2262);
        assert_eq!(pb.prime_pi(10000), 1229); // use existing primes
    }

    #[test]
    fn pb_is_prime_test() {
        // test for is_prime
        let pb = NaiveBuffer::new();

        // some mersenne numbers
        assert!(matches!(
            pb.is_prime(&(2u32.pow(19) - 1), None),
            Primality::Yes
        ));
        assert!(matches!(
            pb.is_prime(&(2u32.pow(23) - 1), None),
            Primality::No
        ));
        assert!(matches!(
            pb.is_prime(&(2u128.pow(89) - 1), None),
            Primality::Probable(_)
        ));

        // test against small prime assertion
        for _ in 0..100 {
            let target = random::<u64>();
            assert_eq!(
                !is_prime64(target),
                matches!(
                    pb.is_prime(&target, Some(PrimalityTestConfig::bpsw())),
                    Primality::No
                )
            );
        }

        // test large numbers
        #[cfg(feature = "num-bigint")]
        {
            let large_primes = [
                // https://golang.org/issue/638
                "18699199384836356663",

                "98920366548084643601728869055592650835572950932266967461790948584315647051443",
                "94560208308847015747498523884063394671606671904944666360068158221458669711639",

                // http://primes.utm.edu/lists/small/small3.html
                "449417999055441493994709297093108513015373787049558499205492347871729927573118262811508386655998299074566974373711472560655026288668094291699357843464363003144674940345912431129144354948751003607115263071543163",
                "230975859993204150666423538988557839555560243929065415434980904258310530753006723857139742334640122533598517597674807096648905501653461687601339782814316124971547968912893214002992086353183070342498989426570593",
                "5521712099665906221540423207019333379125265462121169655563495403888449493493629943498064604536961775110765377745550377067893607246020694972959780839151452457728855382113555867743022746090187341871655890805971735385789993",
                "203956878356401977405765866929034577280193993314348263094772646453283062722701277632936616063144088173312372882677123879538709400158306567338328279154499698366071906766440037074217117805690872792848149112022286332144876183376326512083574821647933992961249917319836219304274280243803104015000563790123",
                // ECC primes: http://tools.ietf.org/html/draft-ladd-safecurves-02
                "3618502788666131106986593281521497120414687020801267626233049500247285301239",                                                                                  // Curve1174: 2^251-9
                "57896044618658097711785492504343953926634992332820282019728792003956564819949",                                                                                 // Curve25519: 2^255-19
                "9850501549098619803069760025035903451269934817616361666987073351061430442874302652853566563721228910201656997576599",                                           // E-382: 2^382-105
                "42307582002575910332922579714097346549017899709713998034217522897561970639123926132812109468141778230245837569601494931472367",                                 // Curve41417: 2^414-17
                "6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151", // E-521: 2^521-1
            ];
            for pstr in large_primes {
                assert!(
                    pb.is_prime(&BigUint::from_str(pstr).unwrap(), None)
                        .probably(),
                    "false negative prime {}",
                    pstr
                )
            }

            let large_composites = [
                "21284175091214687912771199898307297748211672914763848041968395774954376176754",
                "6084766654921918907427900243509372380954290099172559290432744450051395395951",
                "84594350493221918389213352992032324280367711247940675652888030554255915464401",
                "82793403787388584738507275144194252681",

                // Arnault, "Rabin-Miller Primality Test: Composite Numbers Which Pass It",
                // Mathematics of Computation, 64(209) (January 1995), pp. 335-361.
                "1195068768795265792518361315725116351898245581", // strong pseudoprime to prime bases 2 through 29
                // strong pseudoprime to all prime bases up to 200
                "8038374574536394912570796143419421081388376882875581458374889175222974273765333652186502336163960045457915042023603208766569966760987284043965408232928738791850869166857328267761771029389697739470167082304286871099974399765441448453411558724506334092790222752962294149842306881685404326457534018329786111298960644845216191652872597534901",
            ];
            for cstr in large_composites {
                assert!(
                    !pb.is_prime(
                        &BigUint::from_str(cstr).unwrap(),
                        Some(PrimalityTestConfig::bpsw())
                    )
                    .probably(),
                    "false positive prime {}",
                    cstr
                )
            }
        }
    }

    #[test]
    fn pb_factors_test() {
        let pb = NaiveBuffer::new();

        #[cfg(feature = "num-bigint")]
        {
            let m131 = BigUint::from(2u8).pow(131) - 1u8; // m131/263 is a large prime
            let fac = pb.factors(m131, None);
            assert!(matches!(fac, Ok(f) if f.len() == 2));
        }
    }
}
