use std::fs::File;
use std::io::{Write, Error};

use num_prime::factor::{pollard_rho, squfof, one_line, SQUFOF_MULTIPLIERS};
use num_prime::RandPrime;
use rand::random;

fn profile_n(n: u128) -> Vec::<(String, usize)> {    
    let k_squfof: Vec<u16> = SQUFOF_MULTIPLIERS.iter().take(10).cloned().collect();
    let k_oneline: Vec<u16> = vec![1, 360, 480];
    const MAXITER: usize = 1 << 20;

    let mut n_stats = Vec::new();
    
    // pollard rho
    n_stats.push(("pollard_rho1".to_string(), pollard_rho(&n, random(), random(), MAXITER).1));
    n_stats.push(("pollard_rho2".to_string(), pollard_rho(&n, random(), random(), MAXITER).1));

    // squfof
    for &k in &k_squfof {
        let key = format!("squfof_k{}", k);
        if let Some(kn) = n.checked_mul(k as u128) {
            let n = squfof(&n, kn, MAXITER).1;
            n_stats.push((key, n));
        } else {
            n_stats.push((key, MAXITER));
        };
    }

    // one line
    for &k in &k_oneline {
        let key = format!("one_line_k{}", k);
        if let Some(kn) = n.checked_mul(k as u128) {
            let n = one_line(&n, kn, MAXITER).1;
            n_stats.push((key, n));
        } else {
            n_stats.push((key, MAXITER));
        };
    }

    n_stats
}

/// This program try various factorization methods, and log down their iterations number into a csv file
fn main() -> Result<(), Error> {
    let mut rng = rand::thread_rng();
    const REPEATS: u32 = 4;

    let mut n_list = Vec::<(u128, f32)>::new(); // n and bits of n
    let mut stats: Vec<Vec<(String, usize)>> = Vec::new();

    for total_bits in 10..80 {
        for _ in 0..REPEATS {
            let p1: u128 = rng.gen_prime(total_bits / 2, None);
            let p2: u128 = rng.gen_prime_exact(total_bits - (128 - p1.leading_zeros()) as usize, None);
            if p1 == p2 {
                continue;
            }
    
            let n = p1 * p2;
            n_list.push((n, (n as f64).log2() as f32));
            println!("Semiprime: {} = {} * {}", n, p1, p2);
            stats.push(profile_n(n));
        }
    }

    // Log into the CSV file
    let mut fout = File::create("profile_stats.csv")?;
    fout.write(b"n,n_bits")?;
    for k in stats[0].iter().map(|(k, _)| k) {
        fout.write(b",")?;
        fout.write(k.as_bytes())?;
    }

    for ((n, bits), n_stats) in n_list.iter().zip(stats) {
        fout.write(b"\n")?;
        fout.write(n.to_string().as_bytes())?;
        fout.write(b",")?;
        fout.write(bits.to_string().as_bytes())?;
        for (_, v) in n_stats {
            fout.write(b",")?;
            fout.write(v.to_string().as_bytes())?;
        }
    }

    Ok(())
}
