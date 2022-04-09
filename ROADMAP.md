# Roadmap for v0.3
- Address #2 (add `factorize` API for unfallible factorization, and deprecate `factors64`)

# Roadmap for v0.next
- Implement factorization for Gaussian integers (and other quadratic integers?)
- Implement a wrapper supporting fast modular arithmetics and use it to speed up is_prime64 and factors64
- Implement SIQS
- Add benchmarks for factorization & primality test (ref: SageMath benchmark tests)
- Euler totient

# Roadmap for v1
- Support modular integer in a unified API
- Stablize API, determine minimal verison for each dependency
- Support no_std

# Undetermined
- Use NumAssign as trait bounds and see if there's prominent performance improvement
- Support `rug` and `ibig` as backend
- Support `rug` & `primal` or `primesieve-sys` as PrimeBuffer backend
- Support async and multi-thread
