use num_integer::{Integer};
use num_traits::{RefNum, NumRef, FromPrimitive};
use crate::traits::ModInt;

pub fn pollard_rho<T: Integer + FromPrimitive + NumRef + Clone>(target: &T, start: T, offset: T) -> Option<T>
where for<'r> &'r T: RefNum<T> + ModInt<&'r T, &'r T, Output = T> {
    let mut a = start.clone();
    let mut b = start;
    // marker for loop detection, i = tortoise, j = hare
    // see https://www.cnblogs.com/812-xiao-wen/p/10544546.html
    let (mut i, mut j) = (1usize, 2usize);
    while i > 0 {
        i += 1;
        a = (&a).mulm(&a, &target).addm(&offset, &target);
        if a == b { return None; }

        let diff = if b > a { &b - &a } else { &a - &b }; // abs_diff
        let d = diff.gcd(target);
        if d > T::one() && &d < target { return Some(d) }

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

// TODO: ECM, Quadratic sieve / Prime field sieve, SQUFOF
