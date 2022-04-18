#[macro_use]
extern crate criterion;
use criterion::{Criterion, SamplingMode};
use num_prime::{nt_funcs, RandPrime};
#[cfg(feature = "num-primes")]
use num_primes::{Generator, Verification};
use primal_check::miller_rabin;
use number_theory::NumberTheory;

pub fn bench_is_prime(c: &mut Criterion) {
    const N0: u64 = 1_000_000;
    const STEP: usize = 101;
    const N1: u64 = 8_000_000_000; // larger than u32
    const N2: u64 = N1 + N0;

    let numbers = || (1..N0)
        .step_by(STEP)
        .chain((N1..N2).step_by(STEP));

    let mut group = c.benchmark_group("is_prime (u64)");

    group.bench_function("num-prime (this crate)", |b| {
        b.iter(|| {
            numbers()
                .filter(|&n| nt_funcs::is_prime64(n))
                .count()
        })
    });
    #[cfg(feature = "num-primes")]
    group.bench_function("num-primes", |b| {
        b.iter(|| {
            numbers()
                .filter(|&n| Verification::is_prime(&n.into()))
                .count()
        })
    });
    group.bench_function("primal-check", |b| {
        b.iter(|| {
            numbers()
                .filter(|&n| miller_rabin(n))
                .count()
        })
    });
    group.bench_function("number-theory", |b| {
        b.iter(|| {
            numbers()
                .filter(|&n| NumberTheory::is_prime(&n))
                .count()
        })
    });
    // TODO(v0.4.1): why is number-theory faster? difference in modular power?
    // we can add a markdown of performance after figuring out this

    group.finish();
}

pub fn bench_factorization(c: &mut Criterion) {
    const N0: u64 = 1_000_000;
    const STEP: usize = 501;
    const N1: u64 = 8_000_000_000; // larger than u32
    const N2: u64 = N1 + N0;

    let numbers = || (1..N0)
        .step_by(STEP)
        .chain((N1..N2).step_by(STEP));
    let mut group = c.benchmark_group("factorize (u64)");

    group.bench_function("num-prime (this crate)", |b| {
        b.iter(|| {
            numbers()
                .filter(|&n| nt_funcs::factorize64(n as u64).len() > 1)
                .count()
        })
    });
    group.bench_function("number-theory", |b| {
        b.iter(|| {
            numbers()
                .filter(|&n| NumberTheory::factor(&n).len() > 1)
                .count()
        })
    });

    group.finish();
}

pub fn bench_prime_gen(c: &mut Criterion) {
    let mut group = c.benchmark_group("prime generation (256 bits)");
    group.sample_size(10).sampling_mode(SamplingMode::Flat);

    let mut gen = rand::thread_rng();
    group.bench_function("num-prime (this crate)", |b| {
        b.iter(|| -> num_bigint::BigUint { gen.gen_prime(256, None) })
    });
    group.bench_function("num-primes", |b| {
        b.iter(|| Generator::new_prime(256))
    });
    group.finish();

    let mut group = c.benchmark_group("safe prime generation (256 bits)");
    group.sample_size(10).sampling_mode(SamplingMode::Flat);

    let mut gen = rand::thread_rng();
    group.bench_function("num-prime (this crate)", |b| {
        b.iter(|| -> num_bigint::BigUint { gen.gen_safe_prime(256) })
    });
    group.bench_function("num-primes", |b| {
        b.iter(|| Generator::safe_prime(256))
    });
    group.finish();
}

criterion_group!(benches, bench_is_prime, bench_factorization, bench_prime_gen);
criterion_main!(benches);
