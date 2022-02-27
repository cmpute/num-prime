# num-prime

This crate provides utilities for prime number related functionalities:
- Primality testing
  - Deterministic primality check of `u64` integers (using a very fast hashing algorithm)
  - Fermat probable prime test
  - Miller-rabin probable prime test
  - (strong/extra strong) Lucas probable prime test
  - Baillie-PSW test
  - Sophie Germain safe prime test
- Primes generation and indexing
  - A naive implementation of the sieve of Eratosthenes
  - Unified API to support other prime generation backends
  - Generate random (safe) primes
  - Find previous/next prime
- Integer factorization
  - Trial division
  - Pollard's rho algorithm
  - Shanks's square forms factorization (SQUFOF)
  - Fast factorization of `u64` integers
- Number theoretic functions
  - Prime Pi function (number of primes under limit) and its bounds
  - Nth prime and its bounds
  - Moebius function

It's based on the `num` creates and most functions are decently optimized with pre-computed tables.
