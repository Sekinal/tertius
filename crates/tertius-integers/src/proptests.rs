//! Property-based tests for arbitrary precision arithmetic.

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use proptest::prelude::*;

    use crate::{Integer, Rational};

    // Strategy for generating small integers
    fn small_int() -> impl Strategy<Value = i64> {
        -1000i64..1000i64
    }

    // Strategy for generating non-zero integers
    fn non_zero_int() -> impl Strategy<Value = i64> {
        prop_oneof![(-1000i64..=-1i64), (1i64..=1000i64)]
    }

    proptest! {
        // Integer ring axioms

        #[test]
        fn integer_add_commutative(a in small_int(), b in small_int()) {
            let a = Integer::new(a);
            let b = Integer::new(b);
            prop_assert_eq!(a.clone() + b.clone(), b.clone() + a.clone());
        }

        #[test]
        fn integer_add_associative(a in small_int(), b in small_int(), c in small_int()) {
            let a = Integer::new(a);
            let b = Integer::new(b);
            let c = Integer::new(c);
            prop_assert_eq!(
                (a.clone() + b.clone()) + c.clone(),
                a.clone() + (b.clone() + c.clone())
            );
        }

        #[test]
        fn integer_mul_commutative(a in small_int(), b in small_int()) {
            let a = Integer::new(a);
            let b = Integer::new(b);
            prop_assert_eq!(a.clone() * b.clone(), b.clone() * a.clone());
        }

        #[test]
        fn integer_mul_associative(a in small_int(), b in small_int(), c in small_int()) {
            let a = Integer::new(a);
            let b = Integer::new(b);
            let c = Integer::new(c);
            prop_assert_eq!(
                (a.clone() * b.clone()) * c.clone(),
                a.clone() * (b.clone() * c.clone())
            );
        }

        #[test]
        fn integer_distributive(a in small_int(), b in small_int(), c in small_int()) {
            let a = Integer::new(a);
            let b = Integer::new(b);
            let c = Integer::new(c);
            prop_assert_eq!(
                a.clone() * (b.clone() + c.clone()),
                a.clone() * b.clone() + a.clone() * c.clone()
            );
        }

        #[test]
        fn integer_add_identity(a in small_int()) {
            let a = Integer::new(a);
            let zero = Integer::new(0);
            prop_assert_eq!(a.clone() + zero.clone(), a.clone());
            prop_assert_eq!(zero + a.clone(), a);
        }

        #[test]
        fn integer_mul_identity(a in small_int()) {
            let a = Integer::new(a);
            let one = Integer::new(1);
            prop_assert_eq!(a.clone() * one.clone(), a.clone());
            prop_assert_eq!(one * a.clone(), a);
        }

        #[test]
        fn integer_additive_inverse(a in small_int()) {
            let a = Integer::new(a);
            let neg_a = -a.clone();
            let zero = Integer::new(0);
            prop_assert_eq!(a + neg_a, zero);
        }

        // GCD properties

        #[test]
        fn gcd_divides_both(a in non_zero_int(), b in non_zero_int()) {
            let a = Integer::new(a);
            let b = Integer::new(b);
            let g = a.gcd(&b);

            // g should divide both a and b
            let rem_a = a.clone() % g.clone();
            let rem_b = b.clone() % g.clone();
            prop_assert!(rem_a.is_zero());
            prop_assert!(rem_b.is_zero());
        }

        #[test]
        fn gcd_commutative(a in non_zero_int(), b in non_zero_int()) {
            let a = Integer::new(a);
            let b = Integer::new(b);
            prop_assert_eq!(a.gcd(&b), b.gcd(&a));
        }

        // Rational field axioms

        #[test]
        fn rational_add_commutative(
            num_a in small_int(),
            den_a in non_zero_int(),
            num_b in small_int(),
            den_b in non_zero_int()
        ) {
            let a = Rational::from_i64(num_a, den_a);
            let b = Rational::from_i64(num_b, den_b);
            prop_assert_eq!(a.clone() + b.clone(), b.clone() + a.clone());
        }

        #[test]
        fn rational_mul_commutative(
            num_a in small_int(),
            den_a in non_zero_int(),
            num_b in small_int(),
            den_b in non_zero_int()
        ) {
            let a = Rational::from_i64(num_a, den_a);
            let b = Rational::from_i64(num_b, den_b);
            prop_assert_eq!(a.clone() * b.clone(), b.clone() * a.clone());
        }

        #[test]
        fn rational_distributive(
            num_a in small_int(),
            den_a in non_zero_int(),
            num_b in small_int(),
            den_b in non_zero_int(),
            num_c in small_int(),
            den_c in non_zero_int()
        ) {
            let a = Rational::from_i64(num_a, den_a);
            let b = Rational::from_i64(num_b, den_b);
            let c = Rational::from_i64(num_c, den_c);
            prop_assert_eq!(
                a.clone() * (b.clone() + c.clone()),
                a.clone() * b.clone() + a.clone() * c.clone()
            );
        }

        #[test]
        fn rational_multiplicative_inverse(
            num in non_zero_int(),
            den in non_zero_int()
        ) {
            use num_traits::One;
            let a = Rational::from_i64(num, den);
            let inv = a.recip();
            let product = a * inv;
            prop_assert!(product.is_one());
        }

        // ModInt properties

        #[test]
        fn modint_add_commutative(a in 0u64..1000u64, b in 0u64..1000u64) {
            use crate::ModInt;
            const P: u64 = 998244353;
            let a = ModInt::<P>::new(a);
            let b = ModInt::<P>::new(b);
            prop_assert_eq!(a + b, b + a);
        }

        #[test]
        fn modint_mul_commutative(a in 0u64..1000u64, b in 0u64..1000u64) {
            use crate::ModInt;
            const P: u64 = 998244353;
            let a = ModInt::<P>::new(a);
            let b = ModInt::<P>::new(b);
            prop_assert_eq!(a * b, b * a);
        }

        #[test]
        fn modint_inverse(a in 1u64..1000u64) {
            use crate::ModInt;
            const P: u64 = 998244353;
            let a = ModInt::<P>::new(a);
            let inv = a.inv().expect("inverse should exist for non-zero mod prime");
            prop_assert_eq!((a * inv).value(), 1);
        }

        #[test]
        fn modint_fermat_little_theorem(a in 1u64..1000u64) {
            use crate::ModInt;
            const P: u64 = 998244353;
            let a = ModInt::<P>::new(a);
            // a^(p-1) = 1 (mod p) for a != 0
            prop_assert_eq!(a.pow(P - 1).value(), 1);
        }
    }
}
