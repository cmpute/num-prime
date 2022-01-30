use num_integer::{Integer, Roots};
use num_traits::Pow;

/// This trait support unified bit testing for (unsigned) integers
pub trait BitTest {
    fn bits(&self) -> usize;
    fn bit(&self, position: usize) -> bool;
}

/// This enum describes the result of a primality check
pub enum Primality {
    Yes,
    No,
    /// carrying the probability of the number being a prime
    Probable(f32),
}

impl Primality {
    /// Check whether the resule indicates that the number is
    /// (very) probably a prime. Return false only on Primaliy::No
    pub fn probably(self) -> bool {
        match self {
            Primality::No => false, _ => true
        }
    }
}

#[derive(Clone, Copy)]
pub struct PrimalityTestConfig {
    pub sprp_trials: usize,        // number of SPRP test, starting from base 2
    pub sprp_random_trials: usize, // number of SPRP test with random base
    pub slprp_test: bool,
    pub eslprp_test: bool,
}

impl PrimalityTestConfig {
    pub fn default() -> Self {
        Self {
            sprp_trials: 2,
            sprp_random_trials: 2,
            slprp_test: false,
            eslprp_test: false,
        }
    }

    /// Create a configuration for Baillie-PSW test (base 2 SPRP test + SLPRP test)
    pub fn bpsw() -> Self {
        Self {
            sprp_trials: 1,
            sprp_random_trials: 0,
            slprp_test: true,
            eslprp_test: false,
        }
    }

    /// Create a configuration for PSW test (base 2 SPRP + Fibonacci test)
    pub fn psw() {
        todo!()
    } // TODO: implement Fibonacci PRP
}

#[derive(Clone, Copy)]
pub struct FactorizationConfig {
    /// config for test if a
    pub prime_test_config: PrimalityTestConfig,

    /// prime limit of trial division, you also need to reserve the buffer if all primes under the limit are to be tested.
    /// None means using all available primes
    pub td_limit: Option<u64>,

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
            td_limit: Some(1 << 14),
            rho_trials: 4,
            brent_trials: 0,
            pm1_trials: 0,
            pp1_trials: 0,
        }
    }
}

/// Extension on num_integer::Roots to support power check on integers
// FIXME: backport to num_integer (see https://github.com/rust-num/num-traits/issues/233)
pub trait ExactRoots: Roots + Pow<u32, Output = Self> + Clone {
    fn nth_root_exact(&self, n: u32) -> Option<Self> {
        let r = self.nth_root(n);
        if &r.clone().pow(n) == self {
            Some(r)
        } else {
            None
        }
    }
    fn sqrt_exact(&self) -> Option<Self> {
        self.nth_root_exact(2)
    }
    fn cbrt_exact(&self) -> Option<Self> {
        self.nth_root_exact(3)
    }
    fn is_nth_power(&self, n: u32) -> bool {
        self.nth_root_exact(n).is_some()
    }
    fn is_square(&self) -> bool {
        self.sqrt_exact().is_some()
    }
    fn is_cubic(&self) -> bool {
        self.cbrt_exact().is_some()
    }
}

/// This trait describes modular arithmetic on a integer
pub trait ModInt<Rhs = Self, Modulus = Self> {
    type Output;

    /// Return (self + rhs) % m
    fn addm(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (self + rhs) % m
    fn subm(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (self * rhs) % m
    fn mulm(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (self ^ exp) % m
    fn powm(self, exp: Rhs, m: Modulus) -> Self::Output;

    /// Return (-self) % m and make sure the result is in range [0,m)
    fn negm(self, m: Modulus) -> Self::Output;

    /// Calculate inverse module (x such that self*x = 1 mod m)
    fn invm(self, m: Modulus) -> Option<Self::Output>
    where
        Self: Sized;

    /// Return the exponent of factor 2 in the number, usually implemented as trailing_zeros()
    /// This is not directly related to modular arithmetics, but used for implementations of them
    fn fac2(self) -> usize;

    /// Calculate Jacobi Symbol (a|n), where a is self
    fn jacobi(self, n: Modulus) -> i8;

    // TODO: Calculate Kronecker Symbol (a|n), where a is self
    // fn kronecker(self, n: Modulus) -> i8;

    // TODO: Modular sqrt
    // fn sqrtm(self, m: Modulus);
}

/// It's recommended to store at least a bunch of small primes in the buffer
/// to make some of the algorithms more efficient
pub trait PrimeBuffer<'a> {
    // TODO: support indexing?

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

/// This trait implements utility functions for primality checks
/// Reference:
/// - <http://ntheory.org/pseudoprimes.html>
pub trait PrimalityUtils: Integer + Clone {
    /// Test if the integer is a (Fermat) probable prime
    fn is_prp(&self, base: Self) -> bool;

    /// Test if the integer is a strong probable prime (based on miller-rabin test)
    fn is_sprp(&self, base: Self) -> bool;

    /// Test if the integer is a Lucas probable prime
    /// If either of p, q is not specified, then we will use Selfridge's Method A to choose p, q
    fn is_lprp(&self, p: Option<usize>, q: Option<isize>) -> bool;

    /// Test if the integer is a strong Lucas probable prime
    /// If either of p, q is not specified, then we will use Selfridge's Method A to choose p, q
    fn is_slprp(&self, p: Option<usize>, q: Option<isize>) -> bool;

    /// Test if the integer is an extra strong Lucas probable prime
    /// If p is not specified, then first p starting from 3 such that Jacobi symbol is -1 will be chosen, which is sometimes refered as "Method C"
    fn is_eslprp(&self, p: Option<usize>) -> bool;

    // TODO: implement ECPP test
    // https://en.wikipedia.org/wiki/Elliptic_curve_primality

    // TODO: implement is_vprp (Lucas-V probable prime test)
    // https://arxiv.org/pdf/2006.14425.pdf
}
