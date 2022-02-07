#[macro_use]
extern crate criterion;
use criterion::Criterion;
use num_prime::nt_funcs;

pub fn bench_is_prime(c: &mut Criterion) {
    const N: usize = 1_000_000;
    const STEP: usize = 101;
    let mut group = c.benchmark_group("is_prime");

    group.bench_function("64bit", |b| {
        b.iter(|| {
            (1..N)
                .step_by(STEP)
                .filter(|&n| nt_funcs::is_prime64(n as u64))
                .count()
        })
    });
    group.bench_function("default config", |b| {
        b.iter(|| {
            (1..N)
                .step_by(STEP)
                .filter(|&n| nt_funcs::is_prime(&(n as u64), None).probably())
                .count()
        })
    });

    group.finish();
}


pub fn bench_factorization(c: &mut Criterion) {
    const N: usize = 1_000_000;
    const STEP: usize = 501;
    let mut group = c.benchmark_group("factors");

    group.bench_function("64bit", |b| {
        b.iter(|| {
            (1..N)
                .step_by(STEP)
                .filter(|&n| nt_funcs::factors64(n as u64).len() > 1)
                .count()
        })
    });
    group.bench_function("default config", |b| {
        b.iter(|| {
            (1..N)
                .step_by(STEP)
                .filter(|&n| nt_funcs::factors(n as u64, None).is_ok())
                .count()
        })
    });

    group.finish();
}

criterion_group!(benches, bench_is_prime, bench_factorization);
criterion_main!(benches);
