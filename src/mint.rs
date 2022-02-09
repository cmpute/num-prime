//! Wrapper of integer to makes it efficient in modular arithmetics but still have the same
//! API of normal integers.

// TODO (v0.3): Implement basic support for these

/// Integer with fast modular arithmetics support, based on MontgomeryInt
pub struct Mint<T> (Either<T, MontgomeryInt<T>>);

/// Integer with fast modular arithmetics support, based on MontgomeryBigint
pub struct BigMint<T> (Either<Rc<T>, MontgomeryInt<T>>);
