use num_prime::factor::{pollard_rho, squfof, one_line, SQUFOF_MULTIPLIERS};
use num_prime::RandPrime;
use rand::random;

fn main() {
    let mut rng = rand::thread_rng();

    // let p1: u64 = rng.gen_prime(40, None);
    // let p2: u64 = rng.gen_prime(60, None);

    let p1: u64 = 3486784447;
    let p2: u64 = 94143178889;
    // let p2: u64 = 3486784709;

    println!("Semiprime: {} * {}", p1, p2);
    let n = p1 as u128 * p2 as u128;

    // let n: u128 = 133717415095455410877609739380293;

    const MAXITER: usize = 2 << 20;
    for k in SQUFOF_MULTIPLIERS {
        if let Some(kn) = n.checked_mul(k as u128) {
            println!("squfof k={} result: {:?}", k, squfof(&n, kn, MAXITER));
        }
    }
    // println!("one_line k=1 result: {:?}", one_line(&n, n, MAXITER));
    // println!("one_line k=480 result: {:?}", one_line(&n, n * 480, MAXITER));
    // println!("pollard rho result: {:?}", pollard_rho(&n, random(), random(), MAXITER));
}
