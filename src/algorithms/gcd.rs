use crate::big_digit::{BITS, BigDigit, DoubleBigDigit};
use crate::bigint::Sign::*;
use crate::bigint::ToBigInt;
use crate::biguint::IntDigits;
use crate::{BigInt, BigUint};
use alloc::borrow::Cow;
use num_traits::{One, Signed, Zero};

/// XGCD sets z to the greatest common divisor of a and b and returns z.
/// If extended is true, XGCD returns their value such that z = a*x + b*y.
///
/// Allow the inputs a and b to be zero or negative to GCD
/// with the following definitions.
///
/// If x or y are not nil, GCD sets their value such that z = a*x + b*y.
/// Regardless of the signs of a and b, z is always >= 0.
/// If a == b == 0, GCD sets z = x = y = 0.
/// If a == 0 and b != 0, GCD sets z = |b|, x = 0, y = sign(b) * 1.
/// If a != 0 and b == 0, GCD sets z = |a|, x = sign(a) * 1, y = 0.
pub fn xgcd(
    a_in: &BigInt,
    b_in: &BigInt,
    extended: bool,
) -> (BigInt, Option<BigInt>, Option<BigInt>) {
    if a_in.is_zero() && b_in.is_zero() {
        if extended {
            return (0.into(), Some(0.into()), Some(0.into()));
        } else {
            return (0.into(), None, None);
        }
    }

    if a_in.is_zero() {
        if extended {
            let y = if b_in.sign() == Minus {
                -BigInt::one()
            } else {
                BigInt::one()
            };
            return (b_in.abs(), Some(0.into()), Some(y));
        } else {
            return (b_in.abs(), None, None);
        }
    }

    if b_in.is_zero() {
        if extended {
            let x = if a_in.sign() == Minus {
                -BigInt::one()
            } else {
                BigInt::one()
            };
            return (a_in.abs(), Some(x), Some(0.into()));
        } else {
            return (a_in.abs(), None, None);
        }
    }
    lehmer_gcd(a_in, b_in, extended)
}

/// Lehmer GCD sets z to the greatest common divisor of a and b,
/// which both must be != 0, and returns z.
/// If x or y are not nil, their values are set such that z = a*x + b*y.
/// See Knuth, The Art of Computer Programming, Vol. 2, Section 4.5.2, Algorithm L.
fn lehmer_gcd(
    a_in: &BigInt,
    b_in: &BigInt,
    extended: bool,
) -> (BigInt, Option<BigInt>, Option<BigInt>) {
    let mut a = a_in.abs();
    let mut b = b_in.abs();

    let mut ua = if extended { Some(BigInt::one()) } else { None };
    let mut ub = if extended { Some(BigInt::zero()) } else { None };

    let mut q: BigInt = 0.into();
    let mut r: BigInt = 0.into();
    let mut s: BigInt = 0.into();
    let mut t: BigInt = 0.into();

    if a < b {
        core::mem::swap(&mut a, &mut b);
        core::mem::swap(&mut ua, &mut ub);
    }

    while b.len() > 1 {
        let (u0, u1, v0, v1, even) = lehmer_simulate(&a, &b);

        if v0 != 0 {
            lehmer_update(
                &mut a, &mut b, &mut q, &mut r, &mut s, &mut t, u0, u1, v0, v1, even,
            );

            if extended {
                lehmer_update(
                    ua.as_mut().unwrap(),
                    ub.as_mut().unwrap(),
                    &mut q,
                    &mut r,
                    &mut s,
                    &mut t,
                    u0,
                    u1,
                    v0,
                    v1,
                    even,
                );
            }
        } else {
            euclid_update(
                &mut a, &mut b, &mut ua, &mut ub, &mut q, &mut r, &mut s, &mut t, extended,
            );
        }
    }

    if b.len() > 0 {
        if a.len() > 1 {
            euclid_update(
                &mut a, &mut b, &mut ua, &mut ub, &mut q, &mut r, &mut s, &mut t, extended,
            );
        }

        if b.len() > 0 {
            let mut a_word = a.digits()[0];
            let mut b_word = b.digits()[0];

            if extended {
                let mut ua_word: BigDigit = 1;
                let mut ub_word: BigDigit = 0;
                let mut va: BigDigit = 0;
                let mut vb: BigDigit = 1;
                let mut even = true;

                while b_word != 0 {
                    let q = a_word / b_word;
                    let r = a_word % b_word;
                    a_word = b_word;
                    b_word = r;

                    let k = ua_word.wrapping_add(q.wrapping_mul(ub_word));
                    ua_word = ub_word;
                    ub_word = k;

                    let k = va.wrapping_add(q.wrapping_mul(vb));
                    va = vb;
                    vb = k;
                    even = !even;
                }

                t.assign_single_digit(ua_word, if even { Plus } else { Minus });
                s.assign_single_digit(va, if even { Minus } else { Plus });

                if let Some(ua) = ua.as_mut() {
                    t *= &*ua;
                    s *= ub.unwrap();

                    *ua = &t + &s;
                }
            } else {
                while b_word != 0 {
                    let quotient = a_word % b_word;
                    a_word = b_word;
                    b_word = quotient;
                }
            }
            a.digits_mut()[0] = a_word;
        }
    }

    a.normalize();

    let neg_a = a_in.sign() == Minus;

    let y = if let Some(ref mut ua) = ua {
        let mut tmp = a_in * &*ua;

        if neg_a {
            tmp = -tmp;
            *ua = -ua.clone();
        }

        tmp = &a - &tmp;
        tmp = &tmp / b_in;

        Some(tmp)
    } else {
        None
    };

    // Ensure GCD is positive
    if a.sign() == Minus {
        a = -a;
    }

    (a, ua, y)
}

/// Uses the Lehmer algorithm for extended GCD on unsigned values.
pub fn extended_gcd(
    a_in: Cow<'_, BigUint>,
    b_in: Cow<'_, BigUint>,
    extended: bool,
) -> (BigInt, Option<BigInt>, Option<BigInt>) {
    if a_in.is_zero() && b_in.is_zero() {
        if extended {
            return (b_in.to_bigint().unwrap(), Some(0.into()), Some(0.into()));
        } else {
            return (b_in.to_bigint().unwrap(), None, None);
        }
    }

    if a_in.is_zero() {
        if extended {
            return (b_in.to_bigint().unwrap(), Some(0.into()), Some(1.into()));
        } else {
            return (b_in.to_bigint().unwrap(), None, None);
        }
    }

    if b_in.is_zero() {
        if extended {
            return (a_in.to_bigint().unwrap(), Some(1.into()), Some(0.into()));
        } else {
            return (a_in.to_bigint().unwrap(), None, None);
        }
    }

    let a_in = a_in.to_bigint().unwrap();
    let b_in = b_in.to_bigint().unwrap();

    let mut a = a_in.clone();
    let mut b = b_in.clone();

    let mut ua = if extended { Some(BigInt::one()) } else { None };
    let mut ub = if extended { Some(BigInt::zero()) } else { None };

    if a < b {
        core::mem::swap(&mut a, &mut b);
        core::mem::swap(&mut ua, &mut ub);
    }

    let mut q: BigInt = 0.into();
    let mut r: BigInt = 0.into();
    let mut s: BigInt = 0.into();
    let mut t: BigInt = 0.into();

    while b.len() > 1 {
        let (u0, u1, v0, v1, even) = lehmer_simulate(&a, &b);

        if v0 != 0 {
            lehmer_update(
                &mut a, &mut b, &mut q, &mut r, &mut s, &mut t, u0, u1, v0, v1, even,
            );

            if extended {
                lehmer_update(
                    ua.as_mut().unwrap(),
                    ub.as_mut().unwrap(),
                    &mut q,
                    &mut r,
                    &mut s,
                    &mut t,
                    u0,
                    u1,
                    v0,
                    v1,
                    even,
                );
            }
        } else {
            euclid_update(
                &mut a, &mut b, &mut ua, &mut ub, &mut q, &mut r, &mut s, &mut t, extended,
            );
        }
    }

    if b.len() > 0 {
        if a.len() > 1 {
            euclid_update(
                &mut a, &mut b, &mut ua, &mut ub, &mut q, &mut r, &mut s, &mut t, extended,
            );
        }

        if b.len() > 0 {
            let mut a_word = a.digits()[0];
            let mut b_word = b.digits()[0];

            if extended {
                let mut ua_word: BigDigit = 1;
                let mut ub_word: BigDigit = 0;
                let mut va: BigDigit = 0;
                let mut vb: BigDigit = 1;
                let mut even = true;

                while b_word != 0 {
                    let q = a_word / b_word;
                    let r = a_word % b_word;
                    a_word = b_word;
                    b_word = r;

                    let k = ua_word.wrapping_add(q.wrapping_mul(ub_word));
                    ua_word = ub_word;
                    ub_word = k;

                    let k = va.wrapping_add(q.wrapping_mul(vb));
                    va = vb;
                    vb = k;
                    even = !even;
                }

                t.assign_single_digit(ua_word, if even { Plus } else { Minus });
                s.assign_single_digit(va, if even { Minus } else { Plus });

                if let Some(ua) = ua.as_mut() {
                    t *= &*ua;
                    s *= ub.unwrap();

                    *ua = &t + &s;
                }
            } else {
                while b_word != 0 {
                    let quotient = a_word % b_word;
                    a_word = b_word;
                    b_word = quotient;
                }
            }
            a.digits_mut()[0] = a_word;
        }
    }

    a.normalize();

    let y = ua.as_ref().map(|ua| (&a - (&a_in * ua)) / &b_in);

    (a, ua, y)
}

#[inline]
fn lehmer_simulate(a: &BigInt, b: &BigInt) -> (BigDigit, BigDigit, BigDigit, BigDigit, bool) {
    let m = b.len();
    let n = a.len();

    let h = a.digits()[n - 1].leading_zeros();

    let a1: BigDigit = a.digits()[n - 1] << h
        | ((a.digits()[n - 2] as DoubleBigDigit) >> (BITS as u32 - h)) as BigDigit;

    let a2: BigDigit = if n == m {
        b.digits()[n - 1] << h
            | ((b.digits()[n - 2] as DoubleBigDigit) >> (BITS as u32 - h)) as BigDigit
    } else if n == m + 1 {
        ((b.digits()[n - 2] as DoubleBigDigit) >> (BITS as u32 - h)) as BigDigit
    } else {
        0
    };

    let mut even = false;

    let mut u0 = 0;
    let mut u1 = 1;
    let mut u2 = 0;

    let mut v0 = 0;
    let mut v1 = 0;
    let mut v2 = 1;

    let mut a1 = a1;
    let mut a2 = a2;

    while a2 >= v2 && a1.wrapping_sub(a2) >= v1 + v2 {
        let q = a1 / a2;
        let r = a1 % a2;

        a1 = a2;
        a2 = r;

        let k = u1 + q * u2;
        u0 = u1;
        u1 = u2;
        u2 = k;

        let k = v1 + q * v2;
        v0 = v1;
        v1 = v2;
        v2 = k;

        even = !even;
    }

    // Suppress unused variable warnings
    let _ = u0;
    let _ = u2;
    let _ = v0;
    let _ = v2;

    (u0, u1, v0, v1, even)
}

#[allow(clippy::too_many_arguments)]
fn lehmer_update(
    a: &mut BigInt,
    b: &mut BigInt,
    q: &mut BigInt,
    r: &mut BigInt,
    s: &mut BigInt,
    t: &mut BigInt,
    u0: BigDigit,
    u1: BigDigit,
    v0: BigDigit,
    v1: BigDigit,
    even: bool,
) {
    t.assign_single_digit(u0, if even { Plus } else { Minus });
    s.assign_single_digit(v0, if even { Minus } else { Plus });

    *t *= &*a;
    *s *= &*b;

    r.assign_single_digit(u1, if even { Minus } else { Plus });
    q.assign_single_digit(v1, if even { Plus } else { Minus });

    *r *= &*a;
    *q *= &*b;

    *a = &*t + &*s;
    *b = &*r + &*q;
}

#[allow(clippy::too_many_arguments)]
fn euclid_update(
    a: &mut BigInt,
    b: &mut BigInt,
    ua: &mut Option<BigInt>,
    ub: &mut Option<BigInt>,
    q: &mut BigInt,
    r: &mut BigInt,
    s: &mut BigInt,
    t: &mut BigInt,
    extended: bool,
) {
    use num_integer::Integer;

    let (q_new, r_new) = a.div_rem(b);
    *q = q_new;
    *r = r_new;

    core::mem::swap(a, b);
    core::mem::swap(b, r);

    if extended
        && let Some(ub) = ub.as_mut()
        && let Some(ua) = ua.as_mut()
    {
        *t = ub.clone();
        *s = &*ub * &*q;
        *ub = &*ua - &*s;
        *ua = t.clone();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::str::FromStr;

    use num_traits::FromPrimitive;

    #[cfg(feature = "rand")]
    use crate::bigrand::RandBigInt;
    #[cfg(feature = "rand")]
    use num_traits::{One, Zero};
    #[cfg(feature = "rand")]
    use rand::SeedableRng;
    #[cfg(feature = "rand")]
    use rand_xorshift::XorShiftRng;

    #[cfg(feature = "rand")]
    fn extended_gcd_euclid(a: Cow<'_, BigUint>, b: Cow<'_, BigUint>) -> (BigInt, BigInt, BigInt) {
        if a.is_zero() && b.is_zero() {
            return (0.into(), 0.into(), 0.into());
        }

        let (mut s, mut old_s) = (BigInt::zero(), BigInt::one());
        let (mut t, mut old_t) = (BigInt::one(), BigInt::zero());
        let (mut r, mut old_r) = (b.to_bigint().unwrap(), a.to_bigint().unwrap());

        while !r.is_zero() {
            let quotient = &old_r / &r;
            old_r -= &quotient * &r;
            core::mem::swap(&mut old_r, &mut r);
            old_s -= &quotient * &s;
            core::mem::swap(&mut old_s, &mut s);
            old_t -= quotient * &t;
            core::mem::swap(&mut old_t, &mut t);
        }

        (old_r, old_s, old_t)
    }

    #[test]
    fn test_extended_gcd_example() {
        let a = BigUint::from_u32(240).unwrap();
        let b = BigUint::from_u32(46).unwrap();
        let (q, s_k, t_k) = extended_gcd(Cow::Owned(a), Cow::Owned(b), true);

        assert_eq!(q, BigInt::from_i32(2).unwrap());
        assert_eq!(s_k.unwrap(), BigInt::from_i32(-9).unwrap());
        assert_eq!(t_k.unwrap(), BigInt::from_i32(47).unwrap());
    }

    #[test]
    fn test_extended_gcd_example_not_extended() {
        let a = BigUint::from_u32(240).unwrap();
        let b = BigUint::from_u32(46).unwrap();
        let (q, s_k, t_k) = extended_gcd(Cow::Owned(a), Cow::Owned(b), false);

        assert_eq!(q, BigInt::from_i32(2).unwrap());
        assert_eq!(s_k, None);
        assert_eq!(t_k, None);
    }

    #[test]
    fn test_extended_gcd_example_wolfram() {
        let a = BigInt::from_str("-565721958").unwrap();
        let b = BigInt::from_str("4486780496").unwrap();

        let (q, _s_k, _t_k) = xgcd(&a, &b, true);

        assert_eq!(q, BigInt::from(2));
        assert_eq!(_s_k, Some(BigInt::from(-1_090_996_795i64)));
        assert_eq!(_t_k, Some(BigInt::from(-137_559_848i64)));
    }

    #[test]
    fn test_golang_bignum_negative() {
        let gcd_test_cases = [
            ["0", "0", "0", "0", "0"],
            ["7", "0", "1", "0", "7"],
            ["7", "0", "-1", "0", "-7"],
            ["11", "1", "0", "11", "0"],
            ["7", "-1", "-2", "-77", "35"],
            ["935", "-3", "8", "64515", "24310"],
            ["935", "-3", "-8", "64515", "-24310"],
            ["935", "3", "-8", "-64515", "-24310"],
            ["1", "-9", "47", "120", "23"],
            ["7", "1", "-2", "77", "35"],
            ["935", "-3", "8", "64515", "24310"],
            [
                "935000000000000000",
                "-3",
                "8",
                "64515000000000000000",
                "24310000000000000000",
            ],
            [
                "1",
                "-221",
                "22059940471369027483332068679400581064239780177629666810348940098015901108344",
                "98920366548084643601728869055592650835572950932266967461790948584315647051443",
                "991",
            ],
        ];

        for case in &gcd_test_cases {
            let d_case = BigInt::from_str(case[0]).unwrap();
            let x_case = BigInt::from_str(case[1]).unwrap();
            let y_case = BigInt::from_str(case[2]).unwrap();
            let a_case = BigInt::from_str(case[3]).unwrap();
            let b_case = BigInt::from_str(case[4]).unwrap();

            let (_d, _x, _y) = xgcd(&a_case, &b_case, false);
            assert_eq!(_d, d_case);
            assert_eq!(_x, None);
            assert_eq!(_y, None);

            let (_d, _x, _y) = xgcd(&a_case, &b_case, true);
            assert_eq!(_d, d_case);
            assert_eq!(_x.unwrap(), x_case);
            assert_eq!(_y.unwrap(), y_case);
        }
    }

    #[test]
    #[cfg(feature = "rand")]
    fn test_extended_gcd_assumptions() {
        let mut rng = XorShiftRng::from_seed([1u8; 16]);

        for i in 1usize..100 {
            for j in &[1usize, 64, 128] {
                let a = rng.gen_biguint((i * j) as u64);
                let b = rng.gen_biguint((i * j) as u64);
                let (q, s_k, t_k) = extended_gcd(Cow::Borrowed(&a), Cow::Borrowed(&b), true);

                let lhs = BigInt::from_biguint(Plus, a) * &s_k.unwrap();
                let rhs = BigInt::from_biguint(Plus, b) * &t_k.unwrap();

                assert_eq!(q.clone(), &lhs + &rhs, "{} = {} + {}", q, lhs, rhs);
            }
        }
    }

    #[test]
    #[cfg(feature = "rand")]
    fn test_gcd_lehmer_euclid_extended() {
        let mut rng = XorShiftRng::from_seed([1u8; 16]);

        for i in 1usize..80 {
            for j in &[1usize, 16, 24, 64, 128] {
                let a = rng.gen_biguint((i * j) as u64);
                let b = rng.gen_biguint((i * j) as u64);
                let (q, s_k, t_k) = extended_gcd(Cow::Borrowed(&a), Cow::Borrowed(&b), true);

                let expected = extended_gcd_euclid(Cow::Borrowed(&a), Cow::Borrowed(&b));
                assert_eq!(q, expected.0);
                assert_eq!(s_k.unwrap(), expected.1);
                assert_eq!(t_k.unwrap(), expected.2);
            }
        }
    }

    #[test]
    #[cfg(feature = "rand")]
    fn test_gcd_lehmer_euclid_not_extended() {
        let mut rng = XorShiftRng::from_seed([1u8; 16]);

        for i in 1usize..80 {
            for j in &[1usize, 16, 24, 64, 128] {
                let a = rng.gen_biguint((i * j) as u64);
                let b = rng.gen_biguint((i * j) as u64);
                let (q, s_k, t_k) = extended_gcd(Cow::Borrowed(&a), Cow::Borrowed(&b), false);

                let expected = extended_gcd_euclid(Cow::Borrowed(&a), Cow::Borrowed(&b));
                assert_eq!(
                    q, expected.0,
                    "gcd({}, {}) = {} != {}",
                    &a, &b, &q, expected.0
                );
                assert_eq!(s_k, None);
                assert_eq!(t_k, None);
            }
        }
    }
}
