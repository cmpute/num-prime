pub use num_integer::Integer;
use num_traits::{RefNum, NumOps, NumRef, FromPrimitive};
use rand::{Rng, thread_rng};
use rand::distributions::Uniform;
use rand::distributions::uniform::SampleUniform;

/// This trait describes modular arithmetic on a integer
pub trait ModInt<Rhs = Self, Modulus = Self> {    
    type Output;

    /// Return (self * rhs) % m
    fn mul_mod(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (self ^ exp) % m
    fn pow_mod(self, exp: Rhs, m: Modulus) -> Self::Output;
}

// wrapping simple functions
pub trait ArithmeticHelpers {
    fn trailing_zeros(&self) -> usize;
}

/// This trait describes arithmetic functions on a integer
pub trait Arithmetic : Integer + NumOps {
    /// Test if the integer is a strong probable prime
    fn is_sprp(&self, witness: Self) -> bool;

    /// Generate a factor of the integer using Pollard's Rho algorithm
    fn pollard_rho(&self, offset: Self, trials: u32) -> Option<Self>;
}

impl<T: Integer + ArithmeticHelpers + FromPrimitive + NumRef + SampleUniform + Clone> Arithmetic for T
where for<'r> &'r T: RefNum<T> + std::ops::Shr<usize, Output = T> + ModInt<&'r T, &'r T, Output = T>
{
    fn is_sprp(&self, witness: T) -> bool {
        // find 2^shift*u + 1 = n
        let tm1 = self - T::one();
        let shift = tm1.trailing_zeros();
        let u = &tm1 >> shift;

        let mut x = witness.pow_mod(&u, self);
        if x == T::one() || x == tm1 { return true }

        for _ in 0..shift {
            x = (&x).mul_mod(&x, self);
            if x == tm1 { return true }
        }

        x == T::one()
    }

    fn pollard_rho(&self, offset: Self, trials: u32) -> Option<Self> {
        let mut rng = thread_rng();
        let mut trials = trials;
        'trial_loop: while trials > 0 {
            let mut a = rng.sample(Uniform::new(T::from_u8(2u8).unwrap(), self));
            let mut b = a.clone();
            let mut i = 1; let mut j = 2;
            loop {
                i += 1;
                a = ((&a).mul_mod(&a, &self) + &offset) % self;
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
