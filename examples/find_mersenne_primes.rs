use num_prime::nt_funcs::{is_prime, primes};

/// Find all mersenne primes 2^p-1 where p < 128, return a list of p
fn list_mersenne() -> Vec<u64> {
    primes(128)
        .into_iter()
        .filter(|p| is_prime(&(2u128.pow(*p as u32) - 1), None).probably())
        .collect()
}

fn main() {
    println!("Mersenne primes under 2^128:");
    for p in list_mersenne() {
        println!("2^{} - 1", p);
    }
}
