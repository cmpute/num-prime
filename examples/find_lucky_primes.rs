use num_prime::nt_funcs::is_prime64;

/// Find all lucky numbers under limit
/// Referece: <https://en.wikipedia.org/wiki/Lucky_number>
fn list_lucky_numbers(limit: u64) -> Vec<u64> {
    let mut k = 1;
    let mut numbers: Vec<_> = (1..limit).step_by(2).collect();
    while k < numbers.len() / 2 {
        let step = numbers[k] as usize;
        numbers = numbers.into_iter().enumerate().filter(|(i, _)| (i+1) % step != 0).map(|(_, n)| n).collect();
        k += 1;
    }
    numbers
}

/// Find all lucky primes under limit.
/// Reference: <https://prime-numbers.info/article/lucky-primes>
fn list_lucky_primes(limit: u64) -> Vec<u64> {
    list_lucky_numbers(limit).into_iter().filter(|p| is_prime64(*p)).collect()
}

fn main() {
    println!("Lucky primes under 256: {:?}", list_lucky_primes(256));
}
