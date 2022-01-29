//! Implementations for various factorization algorithms

use crate::traits::ModInt;
use num_integer::{Integer, Roots};
use std::collections::BTreeMap;
use num_traits::{FromPrimitive, NumRef, RefNum};

/// Find factors by trial division. The target is guaranteed fully factored
/// only if bound() * bound() > target. The parameter limit sets the max prime to be tried aside from bound()
/// Return factors, and if the target is fully factored, return Ok(residual), otherwise return Err(residual)
pub fn trial_division<I: Iterator<Item = u64>, T: Integer + Clone + Roots + NumRef + FromPrimitive>(
    primes: I,
    target: T,
    limit: Option<u64>,
) -> (BTreeMap<u64, usize>, Result<T, T>)
where
    for<'r> &'r T: RefNum<T>,
{
    let tsqrt: T = num_integer::sqrt(target.clone()) + T::one();
    let limit = if let Some(l) = limit {
        tsqrt.clone().min(T::from_u64(l).unwrap())
    } else {
        tsqrt.clone()
    };

    let mut residual = target;
    let mut result = BTreeMap::new();
    let mut factored = false;
    for (p, pt) in primes.map(|p| (p, T::from_u64(p).unwrap())) {
        if &pt > &tsqrt {
            factored = true;
        }
        if &pt > &limit {
            break;
        }

        while residual.is_multiple_of(&pt) {
            residual = residual / &pt;
            *result.entry(p).or_insert(0) += 1;
        }
        if residual.is_one() {
            factored = true;
            break;
        }
    }

    if factored {
        (result, Ok(residual))
    } else {
        (result, Err(residual))
    }
}


pub fn pollard_rho<T: Integer + FromPrimitive + NumRef + Clone>(
    target: &T,
    start: T,
    offset: T,
) -> Option<T>
where
    for<'r> &'r T: RefNum<T> + ModInt<&'r T, &'r T, Output = T>,
{
    let mut a = start.clone();
    let mut b = start;
    // marker for loop detection, i = tortoise, j = hare
    // see https://www.cnblogs.com/812-xiao-wen/p/10544546.html
    let (mut i, mut j) = (1usize, 2usize);
    while i > 0 {
        i += 1;
        a = (&a).mulm(&a, &target).addm(&offset, &target);
        if a == b {
            return None;
        }

        let diff = if b > a { &b - &a } else { &a - &b }; // abs_diff
        let d = diff.gcd(target);
        if d > T::one() && &d < target {
            return Some(d);
        }

        // when a catches up with b
        if i == j {
            b = a.clone();
            j <<= 1;
        }
    }
    None
}

pub fn pollard_brent() {}
pub fn pollard_pp1() {}
pub fn williams_pp1() {}

// TODO: ECM, Quadratic sieve / Prime field sieve, SQUFOF, Fermat(https://en.wikipedia.org/wiki/Fermat%27s_factorization_method)
