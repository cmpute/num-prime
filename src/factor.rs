//! Implementations for various factorization algorithms.
//!
//! Note general prime number field sieve is not planned to be implemented, since it's too complex
//!
//! See <https://web.archive.org/web/20110331180514/https://diamond.boisestate.edu/~liljanab/BOISECRYPTFall09/Jacobsen.pdf>
//! for a detailed comparison between different factorization algorithms


// XXX: make the factorization method resumable?

use crate::traits::ExactRoots;
use num_integer::{Integer, Roots};
use num_modular::{ModularCoreOps, ModularUnaryOps};
use num_traits::{FromPrimitive, NumRef, RefNum, CheckedAdd};
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
/// 
/// The returned values are the factor and the count of passed iterations.
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
    max_iter: usize,
) -> (Option<T>, usize)
where
    for<'r> &'r T: RefNum<T>,
{
    let mut a = start.clone();
    let mut b = start.clone();
    let mut z = T::one() % target; // accumulator for gcd

    // using Brent's loop detection, i = tortoise, j = hare
    let (mut i, mut j) = (0usize, 1usize);

    // backtracing states
    let mut s = start;
    let mut backtrace = false;

    while i < max_iter {
        i += 1;
        a = a.sqm(&target).addm(&offset, &target);
        if a == b {
            return (None, i);
        }

        // FIXME: optimize abs_diff for montgomery form if we are going to use the abs_diff in the std lib
        let diff = if b > a { &b - &a } else { &a - &b }; // abs_diff
        z = z.mulm(&diff, &target);
        if z.is_zero() {
            // the factor is missed by a combined GCD, do backtracing
            if backtrace {
                // ultimately failed
                return (None, i);
            } else {
                backtrace = true;
                a = std::mem::replace(&mut s, T::one()); // s is discarded
                z = T::one() % target; // clear gcd
                continue;
            }
        }

        // here we check gcd every 2^k steps or 128 steps
        // larger batch size leads to large overhead when backtracing.
        // reference: https://www.cnblogs.com/812-xiao-wen/p/10544546.html
        if i == j || i & 127 == 0 || backtrace {
            let d = z.gcd(target);
            if !d.is_one() && &d != target {
                return (Some(d), i);
            }

            // save state
            s = a.clone();
        }

        // when tortoise catches up with hare, let hare jump to the next stop
        if i == j {
            b = a.clone();
            j <<= 1;
        }
    }

    (None, i)
}

/// This function implements Shanks's square forms factorization (SQUFOF).
///
/// The input is usually multiplied by a multiplier, and the multiplied integer should be put in
/// the `mul_target` argument. The multiplier can be choosen from SQUFOF_MULTIPLIERS, or other square-free odd numbers.
/// The returned values are the factor and the count of passed iterations.
/// 
/// The max iteration can be choosed as 2*n^(1/4), based on Theorem 4.22 from [1].
///
/// Reference: Gower, J., & Wagstaff Jr, S. (2008). Square form factorization.
/// In [1] [Mathematics of Computation](https://homes.cerias.purdue.edu/~ssw/gowerthesis804/wthe.pdf)
/// or [2] [his thesis](https://homes.cerias.purdue.edu/~ssw/gowerthesis804/wthe.pdf)
/// The code is from [3] [Rosetta code](https://rosettacode.org/wiki/Square_form_factorization)
pub fn squfof<T: Integer + NumRef + Clone + ExactRoots + std::fmt::Debug>(target: &T, mul_target: T, max_iter: usize) -> (Option<T>, usize)
where
    for<'r> &'r T: RefNum<T>,
{
    assert!(&mul_target.is_multiple_of(&target), "mul_target should be multiples of target");
    let rd = Roots::sqrt(&mul_target); // root of k*N

    /// Reduction operator for binary quadratic forms. It's equivalent to
    /// the one used in the `num-irrational` crate, in a little different form.
    /// 
    /// This function reduces (a, b, c) = (qm1, p, q), updates qm1 and q, returns new p.
    #[inline]
    fn rho<T: Integer + Clone + NumRef> (rd: &T, p: &T, q: &mut T, qm1: &mut T) -> T where
        for<'r> &'r T: RefNum<T>, {
        let b = (rd + p).div_floor(&*q);
        let new_p = &b * &*q - p;
        let new_q = if p > &new_p {
            &*qm1 + b * (p - &new_p)
        } else {
            &*qm1 - b * (&new_p - p)
        };

        *qm1 = std::mem::replace(q, new_q);
        new_p
    }

    // forward loop, search principal cycle
    let (mut p, mut q, mut qm1) = (rd.clone(), &mul_target - &rd * &rd, T::one());
    if q.is_zero() {
        // shortcut for perfect square
        return (Some(rd), 0);
    }

    for i in 1..max_iter {
        p = rho(&rd, &p, &mut q, &mut qm1);
        if i.is_odd() {
            if let Some(rq) = q.sqrt_exact() {
                let b = (&rd - &p) / &rq;
                let mut u = b * &rq + &p;
                let (mut v, mut vm1) = ((&mul_target - &u * &u) / &rq, rq);

                // backward loop, search ambiguous cycle
                loop {
                    let new_u = rho(&rd, &u, &mut v, &mut vm1);
                    if new_u == u {
                        break;
                    } else {
                        u = new_u
                    }
                }

                let d = target.gcd(&u);
                if d > T::one() && &d < target {
                   return (Some(d), i)
                }
            }
        }
    };
    (None, max_iter)
}

/// Good squfof multipliers sorted by efficiency descendingly, from Dana Jacobsen.
/// 
/// Note: square-free odd numbers are suitable as SQUFOF multipliers
pub const SQUFOF_MULTIPLIERS: [u16; 38] = [
    3 * 5 * 7 * 11, 3 * 5 * 7, 3 * 5 * 7 * 11 * 13, 3 * 5 * 7 * 13, 3 * 5 * 7 * 11 * 17, 3 * 5 * 11,
    3 * 5 * 7 * 17, 3 * 5, 3 * 5 * 7 * 11 * 19, 3 * 5 * 11 * 13, 3 * 5 * 7 * 19, 3 * 5 * 7 * 13 * 17,
    3 * 5 * 13, 3 * 7 * 11, 3 * 7, 5 * 7 * 11, 3 * 7 * 13, 5 * 7,
    3 * 5 * 17, 5 * 7 * 13, 3 * 5 * 19, 3 * 11, 3 * 7 * 17, 3,
    3 * 11 * 13, 5 * 11, 3 * 7 * 19, 3 * 13, 5, 5 * 11 * 13,
    5 * 7 * 19, 5 * 13, 7 * 11, 7, 3 * 17, 7 * 13,
    11, 1
];

/// William Hart's one line factorization algorithm for 64 bit integers.
///
/// The number to be factored could be multiplied by a smooth number (coprime to the target)
/// to speed up, put the multiplied number in the `mul_target` argument. A good multiplier given by Hart is 480.
/// `iters` determine the range for iterating the inner multiplier itself. The returned values are the factor
/// and the count of passed iterations.
/// 
/// 
/// The one line factorization algorithm is especially good at factoring semiprimes with form pq,
/// where p = next_prime(c^a+d1), p = next_prime(c^b+d2), a and b are close, and c, d1, d2 are small integers.
///
/// Reference: Hart, W. B. (2012). A one line factoring algorithm. Journal of the Australian Mathematical Society, 92(1), 61-69. doi:10.1017/S1446788712000146
// TODO: add multipliers preset for one_line method?
pub fn one_line<T: Integer + NumRef + FromPrimitive + ExactRoots + CheckedAdd>(target: &T, mul_target: T, max_iter: usize) -> (Option<T>, usize)
where
    for<'r> &'r T: RefNum<T>, {
    assert!(&mul_target.is_multiple_of(&target), "mul_target should be multiples of target");

    let mut ikn = mul_target.clone();
    for i in 1..max_iter {
        let s = ikn.sqrt() + T::one(); // assuming target is not perfect square
        let m = &s * &s - &ikn;
        if let Some(t) = m.sqrt_exact() {
            let g = target.gcd(&(s - t));
            if !g.is_one() && &g != target {
                return (Some(g), i);
            }
        }

        // prevent overflow
        ikn = if let Some(n) = (&ikn).checked_add(&mul_target) {
            n
        } else {
            return (None, i)
        }
    }
    return (None, max_iter);
}

// TODO: ECM, (self initialize) Quadratic sieve, Lehman's Fermat(https://en.wikipedia.org/wiki/Fermat%27s_factorization_method, n_factor_lehman)
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
        assert_eq!(pollard_rho(&8051u16, 2, 1, 100).0, Some(97));
        assert!(matches!(pollard_rho(&8051u16, random(), 1, 100).0, Some(i) if i == 97 || i == 83));
        assert_eq!(pollard_rho(&455459u32, 2, 1, 100).0, Some(743));

        // Mint test
        for _ in 0..10 {
            let target = random::<u16>() | 1;
            let start = random::<u16>() % target;
            let offset = random::<u16>() % target;

            let expect = pollard_rho(&target, start, offset, 65536);
            let mint_result = pollard_rho(
                &Mint::from(target),
                MontgomeryInt::new(start, target).into(),
                MontgomeryInt::new(offset, target).into(),
                65536
            );
            assert_eq!(
                expect.0,
                mint_result.0.map(|v| v.value())
            );
        }
    }

    #[test]
    fn squfof_test() {
        // case from wikipedia
        assert_eq!(squfof(&11111u32, 11111u32, 100).0, Some(41));

        // cases from https://rosettacode.org/wiki/Square_form_factorization
        let cases: Vec<u64> = vec![
            2501,
            12851,
            13289,
            75301,
            120787,
            967009,
            997417,
            7091569,

            5214317,
            20834839,
            23515517,
            33409583,
            44524219,

            13290059,
            223553581,
            2027651281,
            11111111111,
            100895598169,
            60012462237239,
            287129523414791,
            9007199254740931,
            11111111111111111,
            314159265358979323,
            384307168202281507,
            419244183493398773,
            658812288346769681,
            922337203685477563,
            1000000000000000127,
            1152921505680588799,
            1537228672809128917,

            // this case should success at step 276, from https://rosettacode.org/wiki/Talk:Square_form_factorization
            4558849,
        ];
        for n in cases {
            let d = squfof(&n, n, 40000).0
            .or(squfof(&n, 3*n, 40000).0)
            .or(squfof(&n, 5*n, 40000).0)
            .or(squfof(&n, 7*n, 40000).0)
            .or(squfof(&n, 11*n, 40000).0);
            assert!(matches!(d, Some(_)), "{}", n);
        }
    }

    #[test]
    fn one_line_test() {
        assert_eq!(one_line(&11111u32, 11111u32, 100).0, Some(271));
    }
}
