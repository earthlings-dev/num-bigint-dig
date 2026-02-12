//! Randomization of big integers
#![cfg(feature = "rand")]
#![cfg_attr(docsrs, doc(cfg(feature = "rand")))]

use rand::distr::uniform::{Error, SampleBorrow, SampleUniform, UniformSampler};
use rand::prelude::*;

use crate::BigInt;
use crate::BigUint;
use crate::Sign::*;

use crate::biguint::biguint_from_vec;

use num_integer::Integer;
use num_traits::{ToPrimitive, Zero};

/// A trait for sampling random big integers.
///
/// The `rand` feature must be enabled to use this. See crate-level documentation for details.
pub trait RandBigInt {
    /// Generate a random [`BigUint`] of the given bit size.
    fn gen_biguint(&mut self, bit_size: u64) -> BigUint;

    /// Generate a random [`BigInt`] of the given bit size.
    fn gen_bigint(&mut self, bit_size: u64) -> BigInt;

    /// Generate a random [`BigUint`] less than the given bound. Fails
    /// when the bound is zero.
    fn gen_biguint_below(&mut self, bound: &BigUint) -> BigUint;

    /// Generate a random [`BigUint`] within the given range. The lower
    /// bound is inclusive; the upper bound is exclusive. Fails when
    /// the upper bound is not greater than the lower bound.
    fn gen_biguint_range(&mut self, lbound: &BigUint, ubound: &BigUint) -> BigUint;

    /// Generate a random [`BigInt`] within the given range. The lower
    /// bound is inclusive; the upper bound is exclusive. Fails when
    /// the upper bound is not greater than the lower bound.
    fn gen_bigint_range(&mut self, lbound: &BigInt, ubound: &BigInt) -> BigInt;
}

fn gen_bits<R: Rng + ?Sized>(rng: &mut R, data: &mut [u32], rem: u64) {
    // `fill` is faster than many `random::<u32>` calls
    rng.fill(data);
    if rem > 0 {
        let last = data.len() - 1;
        data[last] >>= 32 - rem;
    }
}

impl<R: Rng + ?Sized> RandBigInt for R {
    cfg_digit!(
        fn gen_biguint(&mut self, bit_size: u64) -> BigUint {
            let (digits, rem) = bit_size.div_rem(&32);
            let len = (digits + (rem > 0) as u64)
                .to_usize()
                .expect("capacity overflow");
            let mut data = vec![0u32; len];
            gen_bits(self, &mut data, rem);
            biguint_from_vec(data)
        }

        fn gen_biguint(&mut self, bit_size: u64) -> BigUint {
            use core::slice;

            let (digits, rem) = bit_size.div_rem(&32);
            let len = (digits + (rem > 0) as u64)
                .to_usize()
                .expect("capacity overflow");
            let native_digits = Integer::div_ceil(&bit_size, &64);
            let native_len = native_digits.to_usize().expect("capacity overflow");
            let mut data = vec![0u64; native_len];
            unsafe {
                // Generate bits in a `&mut [u32]` slice for value stability
                let ptr = data.as_mut_ptr() as *mut u32;
                debug_assert!(native_len * 2 >= len);
                let data = slice::from_raw_parts_mut(ptr, len);
                gen_bits(self, data, rem);
            }
            #[cfg(target_endian = "big")]
            for digit in &mut data {
                // swap u32 digits into u64 endianness
                *digit = (*digit << 32) | (*digit >> 32);
            }
            biguint_from_vec(data)
        }
    );

    fn gen_bigint(&mut self, bit_size: u64) -> BigInt {
        loop {
            // Generate a random BigUint...
            let biguint = self.gen_biguint(bit_size);
            // ...and then randomly assign it a Sign...
            let sign = if biguint.is_zero() {
                // ...except that if the BigUint is zero, we need to try
                // again with probability 0.5. This is because otherwise,
                // the probability of generating a zero BigInt would be
                // double that of any other number.
                if self.random() {
                    continue;
                } else {
                    NoSign
                }
            } else if self.random() {
                Plus
            } else {
                Minus
            };
            return BigInt::from_biguint(sign, biguint);
        }
    }

    fn gen_biguint_below(&mut self, bound: &BigUint) -> BigUint {
        assert!(!bound.is_zero());
        let bits = bound.bits();
        loop {
            let n = self.gen_biguint(bits);
            if n < *bound {
                return n;
            }
        }
    }

    fn gen_biguint_range(&mut self, lbound: &BigUint, ubound: &BigUint) -> BigUint {
        assert!(*lbound < *ubound);
        if lbound.is_zero() {
            self.gen_biguint_below(ubound)
        } else {
            lbound + self.gen_biguint_below(&(ubound - lbound))
        }
    }

    fn gen_bigint_range(&mut self, lbound: &BigInt, ubound: &BigInt) -> BigInt {
        assert!(*lbound < *ubound);
        if lbound.is_zero() {
            BigInt::from(self.gen_biguint_below(ubound.magnitude()))
        } else if ubound.is_zero() {
            lbound + BigInt::from(self.gen_biguint_below(lbound.magnitude()))
        } else {
            let delta = ubound - lbound;
            lbound + BigInt::from(self.gen_biguint_below(delta.magnitude()))
        }
    }
}

/// The back-end implementing rand's [`UniformSampler`] for [`BigUint`].
#[derive(Clone, Debug)]
pub struct UniformBigUint {
    base: BigUint,
    len: BigUint,
}

impl UniformSampler for UniformBigUint {
    type X = BigUint;

    #[inline]
    fn new<B1, B2>(low_b: B1, high_b: B2) -> Result<Self, Error>
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        let low = low_b.borrow();
        let high = high_b.borrow();
        if low >= high {
            return Err(Error::EmptyRange);
        }
        Ok(UniformBigUint {
            len: high - low,
            base: low.clone(),
        })
    }

    #[inline]
    fn new_inclusive<B1, B2>(low_b: B1, high_b: B2) -> Result<Self, Error>
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        let low = low_b.borrow();
        let high = high_b.borrow();
        if low > high {
            return Err(Error::EmptyRange);
        }
        Self::new(low, high + 1u32)
    }

    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Self::X {
        &self.base + rng.gen_biguint_below(&self.len)
    }

    #[inline]
    fn sample_single<R: Rng + ?Sized, B1, B2>(
        low: B1,
        high: B2,
        rng: &mut R,
    ) -> Result<Self::X, Error>
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        let low = low.borrow();
        let high = high.borrow();
        if low >= high {
            return Err(Error::EmptyRange);
        }
        Ok(rng.gen_biguint_range(low, high))
    }
}

impl SampleUniform for BigUint {
    type Sampler = UniformBigUint;
}

/// The back-end implementing rand's [`UniformSampler`] for [`BigInt`].
#[derive(Clone, Debug)]
pub struct UniformBigInt {
    base: BigInt,
    len: BigUint,
}

impl UniformSampler for UniformBigInt {
    type X = BigInt;

    #[inline]
    fn new<B1, B2>(low_b: B1, high_b: B2) -> Result<Self, Error>
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        let low = low_b.borrow();
        let high = high_b.borrow();
        if low >= high {
            return Err(Error::EmptyRange);
        }
        Ok(UniformBigInt {
            len: (high - low).into_parts().1,
            base: low.clone(),
        })
    }

    #[inline]
    fn new_inclusive<B1, B2>(low_b: B1, high_b: B2) -> Result<Self, Error>
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        let low = low_b.borrow();
        let high = high_b.borrow();
        if low > high {
            return Err(Error::EmptyRange);
        }
        Self::new(low, high + 1u32)
    }

    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Self::X {
        &self.base + BigInt::from(rng.gen_biguint_below(&self.len))
    }

    #[inline]
    fn sample_single<R: Rng + ?Sized, B1, B2>(
        low: B1,
        high: B2,
        rng: &mut R,
    ) -> Result<Self::X, Error>
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        let low = low.borrow();
        let high = high.borrow();
        if low >= high {
            return Err(Error::EmptyRange);
        }
        Ok(rng.gen_bigint_range(low, high))
    }
}

impl SampleUniform for BigInt {
    type Sampler = UniformBigInt;
}

/// A random distribution for [`BigUint`] and [`BigInt`] values of a particular bit size.
///
/// The `rand` feature must be enabled to use this. See crate-level documentation for details.
#[derive(Clone, Copy, Debug)]
pub struct RandomBits {
    bits: u64,
}

impl RandomBits {
    #[inline]
    pub fn new(bits: u64) -> RandomBits {
        RandomBits { bits }
    }
}

impl Distribution<BigUint> for RandomBits {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> BigUint {
        rng.gen_biguint(self.bits)
    }
}

impl Distribution<BigInt> for RandomBits {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> BigInt {
        rng.gen_bigint(self.bits)
    }
}

/// A generic trait for generating random primes.
///
/// *Warning*: This is highly dependent on the provided random number generator,
/// to provide actually random primes.
///
/// # Example
#[cfg_attr(feature = "std", doc = " ```")]
#[cfg_attr(not(feature = "std"), doc = " ```ignore")]
/// use num_bigint_dig::RandPrime;
///
/// let mut rng = rand::rng();
/// let p = rng.gen_prime(1024);
/// assert_eq!(p.bits(), 1024);
/// ```
#[cfg(feature = "prime")]
#[cfg_attr(docsrs, doc(cfg(feature = "prime")))]
pub trait RandPrime {
    /// Generate a random prime number with as many bits as given.
    fn gen_prime(&mut self, bits: usize) -> BigUint;
}

/// A list of small, prime numbers that allows us to rapidly
/// exclude some fraction of composite candidates when searching for a random
/// prime. This list is truncated at the point where smallPrimesProduct exceeds
/// a u64. It does not include two because we ensure that the candidates are
/// odd by construction.
#[cfg(feature = "prime")]
const SMALL_PRIMES: [u8; 15] = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];

/// The product of the values in SMALL_PRIMES and allows us
/// to reduce a candidate prime by this number and then determine whether it's
/// coprime to all the elements of SMALL_PRIMES without further BigUint
/// operations.
#[cfg(feature = "prime")]
static SMALL_PRIMES_PRODUCT: std::sync::LazyLock<BigUint> =
    std::sync::LazyLock::new(|| BigUint::from(16_294_579_238_595_022_365u64));

#[cfg(feature = "prime")]
#[cfg_attr(docsrs, doc(cfg(feature = "prime")))]
impl<R: Rng + ?Sized> RandPrime for R {
    fn gen_prime(&mut self, bit_size: usize) -> BigUint {
        use crate::prime::probably_prime;
        use num_traits::ToPrimitive;

        if bit_size < 2 {
            panic!("prime size must be at least 2-bit");
        }

        let mut b = bit_size % 8;
        if b == 0 {
            b = 8;
        }

        let bytes_len = bit_size.div_ceil(8);
        let mut bytes = alloc::vec![0u8; bytes_len];

        loop {
            self.fill_bytes(&mut bytes);
            // Clear bits in the first byte to make sure the candidate has a size <= bits.
            bytes[0] &= ((1u32 << (b as u32)) - 1) as u8;

            // Don't let the value be too small, i.e, set the most significant two bits.
            // Setting the top two bits, rather than just the top bit,
            // means that when two of these values are multiplied together,
            // the result isn't ever one bit short.
            if b >= 2 {
                bytes[0] |= 3u8.wrapping_shl(b as u32 - 2);
            } else {
                // Here b==1, because b cannot be zero.
                bytes[0] |= 1;
                if bytes_len > 1 {
                    bytes[1] |= 0x80;
                }
            }

            // Make the value odd since an even number this large certainly isn't prime.
            bytes[bytes_len - 1] |= 1u8;

            let mut p = BigUint::from_bytes_be(&bytes);
            // must always be a u64, as the SMALL_PRIMES_PRODUCT is a u64
            let rem = (&p % &*SMALL_PRIMES_PRODUCT).to_u64().unwrap();

            'next: for delta in (0u64..1 << 20).step_by(2) {
                let m = rem + delta;

                for prime in &SMALL_PRIMES {
                    if m.is_multiple_of(u64::from(*prime))
                        && (bit_size > 6 || m != u64::from(*prime))
                    {
                        continue 'next;
                    }
                }

                if delta > 0 {
                    p += BigUint::from(delta);
                }

                break;
            }

            // There is a tiny possibility that, by adding delta, we caused
            // the number to be one bit too long. Thus we check bit length here.
            if p.bits() == bit_size as u64 && probably_prime(&p, 20) {
                return p;
            }
        }
    }
}
