# Roadmap for v0.next
- Implement factorization for Gaussian integers (and other quadratic integers?)
- Implement a wrapper supporting fast modular arithmetics and use it to speed up is_prime64 and factors64
- Implement SIQS
- Euler totient

# Roadmap for v1
- Support modular integer in a unified API
- Stablize API, determine minimal verison for each dependency
- Support no_std
- Split the nt_funcs, factor and primality modules into sub files, and put the tables into the related files.

# Undetermined
- Use NumAssign as trait bounds and see if there's prominent performance improvement
- Support `rug` and `ibig` as backend
- Support `rug` & `primal` or `primesieve-sys` as PrimeBuffer backend
- Support async and multi-thread

# Not in plan
- Support number field sieve factorization (this is efficient only for very larg numbers)
