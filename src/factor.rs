//! Implementations for various factorization algorithms
//! 
//! See <https://web.archive.org/web/20110331180514/https://diamond.boisestate.edu/~liljanab/BOISECRYPTFall09/Jacobsen.pdf>
//! for a detailed comparison between different factorization algorithms

use crate::traits::ExactRoots;
use num_integer::{Integer, Roots};
use num_modular::{ModularCoreOps, ModularUnaryOps};
use num_traits::{FromPrimitive, NumRef, RefNum};
use std::collections::BTreeMap;

/// Find factors by trial division, returns a tuple of the found factors and the residual.
///
/// The target is guaranteed fully factored only if bound * bound > target, where bound = max(primes).
/// The parameter limit additionally sets the maximum of primes to be tried.
/// The residual will be Ok(1) or Ok(p) if fully factored.
/// 
/// TODO: implement fast check for small primes with BigInts in the precomputed table, and skip them in this function
pub fn trial_division<
    I: Iterator<Item = u64>,
    T: Integer + Clone + Roots + NumRef + FromPrimitive,
>(
    primes: I,
    target: T,
    limit: Option<u64>,
) -> (BTreeMap<u64, usize>, Result<T, T>)
where
    for<'r> &'r T: RefNum<T>,
{
    let tsqrt: T = Roots::sqrt(&target) + T::one();
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

/// Find factors using Pollard's rho algorithm with Brent's loop detection algorithm
pub fn pollard_rho<
    T: Integer
        + FromPrimitive
        + NumRef
        + Clone
        + for<'r> ModularCoreOps<&'r T, &'r T, Output = T>
        + for<'r> ModularUnaryOps<&'r T, Output = T>,
>(
    target: &T,
    start: T,
    offset: T,
) -> Option<T>
where
    for<'r> &'r T: RefNum<T>,
{
    let mut a = start.clone();
    let mut b = start;
    // using Brent's loop detection, i = tortoise, j = hare
    // XXX: further optimization see https://www.cnblogs.com/812-xiao-wen/p/10544546.html
    let (mut i, mut j) = (1usize, 2usize);
    while i > 0 {
        i += 1;
        a = a.sqm(&target).addm(&offset, &target);
        if a == b {
            return None;
        }

        // FIXME: optimize abs_diff for montgomery form if we are going to use the abs_diff in the std lib
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

/// This function implements Shanks's square forms factorization (SQUFOF). It will assume that target
/// is not a perfect square and the multiplier is square-free.
/// 
/// The multiplier can be choosen from SQUFOF_MULTIPLIERS, or other square-free odd numbers.
/// 
/// Reference: Gower, J., & Wagstaff Jr, S. (2008). Square form factorization.
/// In [Mathematics of Computation](https://homes.cerias.purdue.edu/~ssw/gowerthesis804/wthe.pdf)
/// or [thesis](https://homes.cerias.purdue.edu/~ssw/gowerthesis804/wthe.pdf)
pub fn squfof<T: Integer + NumRef + Clone + ExactRoots>(target: &T, multiplier: T) -> Option<T>
where
    for<'r> &'r T: RefNum<T>,
{
    let kn = multiplier * target;

    // the strategy of limiting iterations is from GNU factor
    let s = Roots::sqrt(&kn);
    let two = T::one() + T::one();
    let max_iter = &two * Roots::sqrt(&(&two * &s));

    // forward
    let p0 = s;
    let mut pm1 = p0.clone();
    let mut p; // to be initialized in the first iteration
    let mut qm1 = T::one();
    let mut q = &kn - &p0 * &p0;
    let mut i = T::one();
    let qsqrt = loop {
        let b = (&p0 + &pm1) / &q;
        p = &b * &q - &pm1;
        let qnext = if pm1 > p {
            &qm1 + &b * (&pm1 - &p)
        } else {
            &qm1 - &b * (&p - &pm1)
        };
        if i.is_odd() {
            if let Some(v) = qnext.sqrt_exact() {
                break v;
            }
        }

        pm1 = p;
        qm1 = q;
        q = qnext;
        i = i + T::one();

        if i == max_iter {
            return None;
        }
    };

    // backward
    let b0 = (&p0 - &p) / &qsqrt;
    pm1 = &b0 * &qsqrt + &p;
    qm1 = qsqrt;
    q = (&kn - &pm1 * &pm1) / &qm1;

    loop {
        let b = (&p0 + &pm1) / &q;
        p = &b * &q - &pm1;
        if p == pm1 {
            break;
        }

        let qnext = if pm1 > p {
            &qm1 + &b * (&pm1 - &p)
        } else {
            &qm1 - &b * (&p - &pm1)
        };
        pm1 = p;
        qm1 = q;
        q = qnext;
    }

    let d = target.gcd(&p);
    if d > T::one() && &d < target {
        Some(d)
    } else {
        None
    }
}

// Square-free even numbers are suitable as SQUFOF multipliers 
pub const SQUFOF_MULTIPLIERS: [u16; 16] = [
    1,
    3,
    5,
    7,
    11,
    3 * 5,
    3 * 7,
    3 * 11,
    5 * 7,
    5 * 11,
    7 * 11,
    3 * 5 * 7,
    3 * 5 * 11,
    3 * 7 * 11,
    5 * 7 * 11,
    3 * 5 * 7 * 11,
];

// TODO(v0.3.3): implement one line factorization and its optimization
// REF:  doi:10.1017/S1446788712000146
//      https://math.mit.edu/research/highschool/primes/materials/2019/Gopalakrishna.pdf

// TODO: ECM, Quadratic sieve / Prime field sieve, Fermat(https://en.wikipedia.org/wiki/Fermat%27s_factorization_method)
// REF: https://pypi.org/project/primefac/
//      http://flintlib.org/doc/ulong_extras.html#factorisation
//      https://github.com/zademn/facto-rs/
//      https://github.com/elmomoilanen/prime-factorization
//      https://cseweb.ucsd.edu/~ethome/teaching/2022-cse-291-14/
fn pollard_pp1() {}
fn williams_pp1() {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mint::Mint;
    use num_modular::MontgomeryInt;
    use rand::random;

    #[test]
    fn pollard_rho_test() {
        assert_eq!(pollard_rho(&8051u16, 2, 1), Some(97));
        assert!(matches!(pollard_rho(&8051u16, random(), 1), Some(i) if i == 97 || i == 83));
        assert_eq!(pollard_rho(&455459u32, 2, 1), Some(743));

        // Mint test
        for _ in 0..10 {
            let target = random::<u16>() | 1;
            let start = random::<u16>() % target;
            let offset = random::<u16>() % target;
            assert_eq!(pollard_rho(&target, start, offset),
                       pollard_rho(&Mint::from(target),
                       MontgomeryInt::new(start, target).into(),
                       MontgomeryInt::new(offset, target).into()
                ).map(|v| v.value()));
        }
    }

    #[test]
    fn squfof_test() {
        assert!(matches!(squfof(&11111u32, 1), Some(41)));
    }
}
