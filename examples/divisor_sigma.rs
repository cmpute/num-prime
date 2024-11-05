use num_prime::nt_funcs::factorize64;

/// Return all divisors of the target
fn divisors(target: u64) -> Vec<u64> {
    let factors = factorize64(target);
    let mut result = Vec::with_capacity(factors.values().map(|e| e + 1).product());
    result.push(1);

    for (p, e) in factors {
        // the new results contain all previous divisors multiplied by p, p^2, .., p^e
        let mut new_result = Vec::with_capacity(result.len() * e);
        for i in 1..=(e as u32) {
            new_result.extend(result.iter().map(|f| f * p.pow(i)));
        }
        result.append(&mut new_result);
    }
    result
}

/// Calculate the divisor sigma function `σ_z(n`) on the target
/// Reference: <https://en.wikipedia.org/wiki/Divisor_function>
fn divisor_sigma(target: u64, z: u32) -> u64 {
    divisors(target).into_iter().map(|d| d.pow(z)).sum()
}

fn main() {
    println!("Divisor sigma with z=1 of numbers from 10 to 99:");
    for i in 10..100 {
        println!("{}: {:?}", i, divisor_sigma(i, 1));
    }
}
