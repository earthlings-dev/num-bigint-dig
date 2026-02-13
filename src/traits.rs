use alloc::borrow::Cow;

use num_integer::Integer;
use num_traits::Signed;

use crate::algorithms::{extended_gcd, mod_inverse, xgcd};
use crate::{BigInt, BigUint};

/// Generic trait for modular multiplicative inverse.
///
/// Computes the [modular multiplicative inverse](https://en.wikipedia.org/wiki/Modular_multiplicative_inverse)
/// of an integer *a* modulo *m*.
///
/// Returns `None` if the inverse does not exist (i.e., `gcd(a, m) != 1`).
pub trait ModInverse<R: Sized>: Sized {
    /// The output type of the modular inverse.
    type Output: Sized;

    /// Returns the modular inverse of `self` modulo `m`, or `None` if it does not exist.
    fn mod_inverse(self, m: R) -> Option<Self::Output>;
}

/// Generic trait for the extended Euclidean algorithm.
///
/// Computes the [extended GCD](https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm),
/// returning `(gcd, x, y)` such that `self * x + other * y = gcd`.
pub trait ExtendedGcd<R: Sized>: Sized {
    /// Returns `(gcd, x, y)` such that `self * x + other * y = gcd`.
    fn extended_gcd(self, other: R) -> (BigInt, BigInt, BigInt);
}

// --- ModInverse impls ---

impl ModInverse<&BigUint> for BigUint {
    type Output = BigInt;

    fn mod_inverse(self, m: &BigUint) -> Option<BigInt> {
        mod_inverse(Cow::Owned(self), Cow::Borrowed(m))
    }
}

impl ModInverse<BigUint> for BigUint {
    type Output = BigInt;

    fn mod_inverse(self, m: BigUint) -> Option<BigInt> {
        mod_inverse(Cow::Owned(self), Cow::Owned(m))
    }
}

impl ModInverse<&BigUint> for BigInt {
    type Output = BigInt;

    fn mod_inverse(self, m: &BigUint) -> Option<BigInt> {
        if self.is_negative() {
            let v = self
                .mod_floor(&BigInt::from(m.clone()))
                .to_biguint()
                .unwrap();
            mod_inverse(Cow::Owned(v), Cow::Borrowed(m))
        } else {
            mod_inverse(Cow::Owned(self.into_parts().1), Cow::Borrowed(m))
        }
    }
}

impl ModInverse<&BigInt> for BigInt {
    type Output = BigInt;

    fn mod_inverse(self, m: &BigInt) -> Option<BigInt> {
        let modulus = m.magnitude().clone();

        if self.is_negative() {
            let v = self
                .mod_floor(&BigInt::from(modulus.clone()))
                .to_biguint()
                .unwrap();
            mod_inverse(Cow::Owned(v), Cow::Owned(modulus))
        } else {
            mod_inverse(Cow::Owned(self.into_parts().1), Cow::Owned(modulus))
        }
    }
}

impl ModInverse<BigInt> for BigInt {
    type Output = BigInt;

    fn mod_inverse(self, m: BigInt) -> Option<BigInt> {
        if self.is_negative() {
            let v = self.mod_floor(&m).to_biguint().unwrap();
            mod_inverse(Cow::Owned(v), Cow::Owned(m.into_parts().1))
        } else {
            mod_inverse(
                Cow::Owned(self.into_parts().1),
                Cow::Owned(m.into_parts().1),
            )
        }
    }
}

// --- ExtendedGcd impls ---

impl ExtendedGcd<&BigUint> for BigUint {
    fn extended_gcd(self, other: &BigUint) -> (BigInt, BigInt, BigInt) {
        let (gcd, x, y) = extended_gcd(Cow::Owned(self), Cow::Borrowed(other), true);
        (gcd, x.unwrap(), y.unwrap())
    }
}

impl ExtendedGcd<BigUint> for BigUint {
    fn extended_gcd(self, other: BigUint) -> (BigInt, BigInt, BigInt) {
        let (gcd, x, y) = extended_gcd(Cow::Owned(self), Cow::Owned(other), true);
        (gcd, x.unwrap(), y.unwrap())
    }
}

impl ExtendedGcd<&BigInt> for BigInt {
    fn extended_gcd(self, other: &BigInt) -> (BigInt, BigInt, BigInt) {
        let (gcd, x, y) = xgcd(&self, other, true);
        (gcd, x.unwrap(), y.unwrap())
    }
}

impl ExtendedGcd<BigInt> for BigInt {
    fn extended_gcd(self, other: BigInt) -> (BigInt, BigInt, BigInt) {
        let (gcd, x, y) = xgcd(&self, &other, true);
        (gcd, x.unwrap(), y.unwrap())
    }
}
