pub use num_integer::{Integer, Roots};
use num_traits::{RefNum, NumOps, NumRef, FromPrimitive};
use rand::{Rng, thread_rng};
use rand::distributions::Uniform;
use rand::distributions::uniform::SampleUniform;

/// Extension on num_integer::Roots to support power check on integers
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

    /// Return the exponent of factor 2 in the number, usually implemented as trailing_zeros()
    fn fac2(self) -> usize;
    
    /// Calculate Jacobi Symbol (a|n), where a is self
    fn jacobi(self, n: Modulus) -> i8;

    /// Calculate inverse module (x such that self*x = 1 mod m)
    fn invm(self, m: Modulus) -> Option<Self::Output> where Self: Sized;
}

/// This trait describes number theoretic functions on a integer
/// Reference:
/// - <http://ntheory.org/pseudoprimes.html>
/// - <http://www.numbertheory.org/gnubc/bc_programs.html>
pub trait NumberTheoretic : Integer + NumOps + FromPrimitive + NumRef + SampleUniform + Clone {
    /// Test if the integer is a (Fermat) probable prime
    fn is_prp(&self, base: Self) -> bool;

    /// Test if the integer is a strong probable prime (based on miller-rabin test)
    fn is_sprp(&self, base: Self) -> bool;

    // TODO: implement is_slprp (Strong Lucas Probable Prime)
    // https://en.wikipedia.org/wiki/Lucas_pseudoprime

    // TODO: implement ECPP test?
    // https://en.wikipedia.org/wiki/Elliptic_curve_primality

    /// Generate a factor of the integer using Pollard's Rho algorithm
    /// TODO(v0.0.4): remove dependency on SampleUniform?
    fn pollard_rho(&self, offset: Self, trials: usize) -> Option<Self>;
}

// TODO: implement other utilities in https://gmplib.org/manual/Number-Theoretic-Functions

impl<T: Integer + FromPrimitive + NumRef + SampleUniform + Clone> NumberTheoretic for T
where for<'r> &'r T: RefNum<T> + std::ops::Shr<usize, Output = T> + ModInt<&'r T, &'r T, Output = T>
{
    fn is_prp(&self, base: Self) -> bool {
        let tm1 = self - Self::one();
        base.powm(&tm1, self).is_one()
    }

    fn is_sprp(&self, base: T) -> bool {
        // find 2^shift*u + 1 = n
        let tm1 = self - T::one();
        let shift = tm1.fac2();
        let u = &tm1 >> shift;

        let mut x = base.powm(&u, self);
        if x == T::one() || x == tm1 { return true }

        for _ in 0..shift {
            x = (&x).mulm(&x, self);
            if x == tm1 { return true }
        }

        x == T::one()
    }

    fn pollard_rho(&self, offset: Self, trials: usize) -> Option<Self> {
        let mut rng = thread_rng();
        let mut trials = trials;
        'trial_loop: while trials > 0 {
            let mut a = rng.sample(Uniform::new(T::from_u8(2u8).unwrap(), self));
            let mut b = a.clone();
            let mut i = 1; let mut j = 2;
            loop {
                i += 1;
                a = ((&a).mulm(&a, &self) + &offset) % self;
                if a == b {
                    trials -= 1;
                    continue 'trial_loop
                }
                let diff = if b > a { &b - &a } else { &a - &b }; // abs_diff
                let d = diff.gcd(self);
                if T::one() < d && d < a {
                    return Some(d)
                }
                if i == j {
                    b = a.clone();
                    j <<= 1;
                }
            }
        }
        None
    }
}
