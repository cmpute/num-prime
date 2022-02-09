# Roadmap for v0.2
- Implement SIQS
- Add benchmarks for factorization & primality test (ref: SageMath benchmark tests)
- Support no_std
- Use NumAssign as trait bounds and see if there's prominent performance improvement

# Roadmap for v1
- Support `rug` and `ibig` as backend
- Stablize API
- Euler totient
- Implement a wrapper supporting fast modular arithmetics. (Maybe implement it as `Eiter<Integer, ModularInteger>`)
- [?] Support rug & primal or primesieve-sys as backend
- [?] Support async and multi-thread
