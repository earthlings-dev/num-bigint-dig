//! Low-level algorithms for big integer arithmetic.
//!
//! This module re-exports the internal arithmetic primitives used by
//! [`BigUint`](crate::BigUint) and [`BigInt`](crate::BigInt), making them
//! available for direct use on digit slices.

#![allow(clippy::many_single_char_names)]

// Re-export arithmetic primitives from their canonical locations.
//
// The implementations live inside the `biguint` submodules where they are
// used by the operator trait impls. We simply widen their visibility here.

// --- addition ---
pub use crate::biguint::addition::{__add2, adc, add2};

// --- subtraction ---
pub use crate::biguint::subtraction::{__sub2rev, sbb, sub2, sub2rev};

// --- multiplication ---
pub use crate::biguint::multiplication::{
    mac_digit, mac_with_carry, mac3, mul3, scalar_mul, sub_sign,
};

// --- division ---
pub use crate::biguint::division::{div_rem, div_rem_digit, div_rem_ref};

// --- shift ---
pub use crate::biguint::shift::{biguint_shl, biguint_shr};

// --- comparison ---
pub use crate::biguint::cmp_slice;

// --- montgomery modular exponentiation ---
pub use crate::biguint::monty::monty_modpow;

// --- crypto algorithms ---
mod bits;
mod gcd;
mod jacobi;
mod mod_inverse;

pub use self::bits::*;
pub use self::gcd::*;
pub use self::jacobi::*;
pub use self::mod_inverse::*;
