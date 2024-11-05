use num_prime::nt_funcs::factorize64;

/// Calculate the (small) prime omega function ω(n) on the target
/// Reference: <https://en.wikipedia.org/wiki/Prime_omega_function>
fn prime_omega(target: u64) -> usize {
    factorize64(target).len()
}

/// Calculate the (big) prime omega function Ω(n) on the target
/// Reference: <https://en.wikipedia.org/wiki/Prime_omega_function>
#[allow(non_snake_case)]
fn prime_Omega(target: u64) -> usize {
    factorize64(target).into_values().sum()
}

fn main() {
    println!("Prime omega of numbers from 10 to 99:");
    for i in 10..100 {
        println!("{}: ω={}, Ω={}", i, prime_omega(i), prime_Omega(i));
    }
}
