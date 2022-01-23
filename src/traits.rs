pub use num_integer::Integer;
use num_traits::{RefNum, NumOps, NumRef, FromPrimitive};
use rand::{Rng, thread_rng};
use rand::distributions::Uniform;
use rand::distributions::uniform::SampleUniform;

/// This trait describes modular arithmetic on a integer
pub trait ModInt<Rhs = Self, Modulus = Self> {    
    type Output;

    /// Return (self * rhs) % m
    fn mulm(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (self ^ exp) % m
    fn powm(self, exp: Rhs, m: Modulus) -> Self::Output;

    /// Return the exponent of factor 2 in the number, usually implemented as trailing_zeros()
    fn fac2(self) -> usize;
}

/// This trait describes arithmetic functions on a integer
pub trait PrimeArithmetic : Integer + NumOps {
    /// Test if the integer is a strong (Fermat) probable prime
    fn is_sprp(&self, witness: Self) -> bool;

    // TODO: implement is_slprp (Strong Lucas Probable Prime)
    // https://en.wikipedia.org/wiki/Lucas_pseudoprime

    // TODO: implement ECPP test?
    // https://en.wikipedia.org/wiki/Elliptic_curve_primality

    /// Generate a factor of the integer using Pollard's Rho algorithm
    fn pollard_rho(&self, offset: Self, trials: usize) -> Option<Self>;
}

// TODO: implement other utilities in https://gmplib.org/manual/Number-Theoretic-Functions
// TODO: add invm (https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/, https://github.com/JohnHammond/primefac_fork)

impl<T: Integer + FromPrimitive + NumRef + SampleUniform + Clone> PrimeArithmetic for T
where for<'r> &'r T: RefNum<T> + std::ops::Shr<usize, Output = T> + ModInt<&'r T, &'r T, Output = T>
{
    // TODO: implement these functions as default implementation, so that other types can implement themselves
    fn is_sprp(&self, witness: T) -> bool {
        // find 2^shift*u + 1 = n
        let tm1 = self - T::one();
        let shift = tm1.fac2();
        let u = &tm1 >> shift;

        let mut x = witness.powm(&u, self);
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
