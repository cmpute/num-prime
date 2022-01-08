// XXX: implement streaming prime sieve, like `primal` crate

mod buffer;
mod traits;
mod integer;

pub use buffer::PrimeBuffer;
pub use traits::{Arithmetic, ModInt};
