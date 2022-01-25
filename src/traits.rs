use num_integer::{Integer, Roots};

/// Extension on num_integer::Roots to support power check on integers
// TODO: backport to num_integer
pub trait ExactRoot : Roots {
    fn is_nth_power(&self) -> bool;
    fn is_square(&self) -> bool;
    fn is_cubic(&self) -> bool;
    fn nth_root_exact(&self) -> Option<Self>;
    fn sqrt_exact(&self) -> Option<Self>;
    fn cbrt_exact(&self) -> Option<Self>;
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

    /// Calculate inverse module (x such that self*x = 1 mod m)
    fn invm(self, m: Modulus) -> Option<Self::Output> where Self: Sized;

    /// Return the exponent of factor 2 in the number, usually implemented as trailing_zeros()
    /// This is not directly related to modular arithmetics, but used for implementations of them
    fn fac2(self) -> usize;
    
    /// Calculate Jacobi Symbol (a|n), where a is self
    fn jacobi(self, n: Modulus) -> i8;

    // TODO: Calculate Kronecker Symbol (a|n), where a is self
    // fn kronecker(self, n: Modulus) -> i8;
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

/// This trait implements utility functions for primality check and factorization
/// Reference:
/// - <http://ntheory.org/pseudoprimes.html>
/// - <http://www.numbertheory.org/gnubc/bc_programs.html>
pub trait PrimalityUtils : Integer + Clone {
    /// Test if the integer is a (Fermat) probable prime
    fn is_prp(&self, base: Self) -> bool;

    /// Test if the integer is a strong probable prime (based on miller-rabin test)
    fn is_sprp(&self, base: Self) -> bool;

    /// Test if the integer is a strong Lucas probable prime
    fn is_slprp(&self, p: Self, q: Self) -> bool;

    /// Test if the integer is a strong Lucas probable prime with P, Q defined by Selfridge's Method
    fn is_sslprp(&self) -> bool;
    
    // TODO: implement ECPP test?
    // https://en.wikipedia.org/wiki/Elliptic_curve_primality
}
