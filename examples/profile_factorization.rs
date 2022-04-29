use num_prime::factor::{pollard_rho, squfof, one_line};
use num_prime::RandPrime;
use rand::random;

fn main() {
    let mut rng = rand::thread_rng();


    let p1: u64 = rng.gen_prime(40, None);
    let p2: u64 = rng.gen_prime(60, None);

    // let p1: u64 = 3486784447;
    // let p2: u64 = 94143178889;
    // let p2: u64 = 3486784709;

    let n = p1 as u128 * p2 as u128;
    println!("Semiprime: {} * {}", p1, p2);
    println!("pollard rho result: {:?}", pollard_rho(&n, random(), random(), 400000));
    println!("squfof k=1 result: {:?}", squfof(&n, n, 40000));
    println!("squfof k=3*5*7*11 result: {:?}", squfof(&n, n * 3 * 5 * 7 * 11, 40000));
    println!("one_line k=1 result: {:?}", one_line(&n, n, 400000));
    println!("one_line k=480 result: {:?}", one_line(&n, n * 480, 400000));
}
