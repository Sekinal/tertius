//! Property-based tests for polynomial arithmetic.

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use crate::dense::DensePoly;
    use tertius_rings::rationals::Q;

    // Strategy for generating small rational coefficients
    fn small_coeff() -> impl Strategy<Value = Q> {
        (-100i64..100i64).prop_map(Q::from_integer)
    }

    // Strategy for generating small polynomials (degree 0-4)
    fn small_poly() -> impl Strategy<Value = DensePoly<Q>> {
        proptest::collection::vec(small_coeff(), 1..=5).prop_map(DensePoly::new)
    }

    // Strategy for generating non-zero polynomials
    fn nonzero_poly() -> impl Strategy<Value = DensePoly<Q>> {
        small_poly().prop_filter("polynomial must be non-zero", |p| !p.is_zero())
    }

    proptest! {
        // Polynomial ring axioms

        #[test]
        fn poly_add_commutative(a in small_poly(), b in small_poly()) {
            prop_assert_eq!(a.add(&b), b.add(&a));
        }

        #[test]
        fn poly_add_associative(a in small_poly(), b in small_poly(), c in small_poly()) {
            prop_assert_eq!(a.add(&b).add(&c), a.add(&b.add(&c)));
        }

        #[test]
        fn poly_mul_commutative(a in small_poly(), b in small_poly()) {
            prop_assert_eq!(a.mul(&b), b.mul(&a));
        }

        #[test]
        fn poly_mul_associative(a in small_poly(), b in small_poly(), c in small_poly()) {
            prop_assert_eq!(a.mul(&b).mul(&c), a.mul(&b.mul(&c)));
        }

        #[test]
        fn poly_distributive(a in small_poly(), b in small_poly(), c in small_poly()) {
            // a * (b + c) = a * b + a * c
            let left = a.mul(&b.add(&c));
            let right = a.mul(&b).add(&a.mul(&c));
            prop_assert_eq!(left, right);
        }

        #[test]
        fn poly_add_identity(a in small_poly()) {
            let zero = DensePoly::zero();
            prop_assert_eq!(a.add(&zero), a.clone());
            prop_assert_eq!(zero.add(&a), a);
        }

        #[test]
        fn poly_mul_identity(a in small_poly()) {
            let one = DensePoly::one();
            prop_assert_eq!(a.mul(&one), a.clone());
            prop_assert_eq!(one.mul(&a), a);
        }

        #[test]
        fn poly_mul_zero(a in small_poly()) {
            let zero = DensePoly::zero();
            prop_assert!(a.mul(&zero).is_zero());
            prop_assert!(zero.mul(&a).is_zero());
        }

        #[test]
        fn poly_additive_inverse(a in small_poly()) {
            let neg_a = a.neg();
            let sum = a.add(&neg_a);
            prop_assert!(sum.is_zero());
        }

        // Degree properties

        #[test]
        fn poly_mul_degree(a in nonzero_poly(), b in nonzero_poly()) {
            // deg(a * b) = deg(a) + deg(b) for non-zero polynomials
            let product = a.mul(&b);
            if !product.is_zero() {
                prop_assert_eq!(product.degree(), a.degree() + b.degree());
            }
        }

        #[test]
        fn poly_add_degree_bound(a in small_poly(), b in small_poly()) {
            // deg(a + b) <= max(deg(a), deg(b))
            let sum = a.add(&b);
            prop_assert!(sum.degree() <= a.degree().max(b.degree()));
        }

        // Evaluation property

        #[test]
        fn poly_eval_add(a in small_poly(), b in small_poly(), x in small_coeff()) {
            // (a + b)(x) = a(x) + b(x)
            let sum = a.add(&b);
            prop_assert_eq!(sum.eval(&x), a.eval(&x) + b.eval(&x));
        }

        #[test]
        fn poly_eval_mul(a in small_poly(), b in small_poly(), x in small_coeff()) {
            // (a * b)(x) = a(x) * b(x)
            let product = a.mul(&b);
            prop_assert_eq!(product.eval(&x), a.eval(&x) * b.eval(&x));
        }

        // Karatsuba vs schoolbook equivalence

        #[test]
        fn karatsuba_matches_schoolbook(
            a_coeffs in proptest::collection::vec(-10i64..10i64, 1..=8),
            b_coeffs in proptest::collection::vec(-10i64..10i64, 1..=8)
        ) {
            let a = DensePoly::new(a_coeffs.into_iter().map(Q::from_integer).collect());
            let b = DensePoly::new(b_coeffs.into_iter().map(Q::from_integer).collect());

            // Use schoolbook for reference
            let schoolbook_result = schoolbook_mul(&a, &b);
            let karatsuba_result = a.mul(&b);

            prop_assert_eq!(karatsuba_result, schoolbook_result);
        }
    }

    // Helper function for reference schoolbook multiplication
    fn schoolbook_mul(a: &DensePoly<Q>, b: &DensePoly<Q>) -> DensePoly<Q> {
        use tertius_rings::traits::Ring;

        if a.is_zero() || b.is_zero() {
            return DensePoly::zero();
        }

        let mut result = vec![Q::zero(); a.coeffs().len() + b.coeffs().len() - 1];
        for (i, ai) in a.coeffs().iter().enumerate() {
            for (j, bj) in b.coeffs().iter().enumerate() {
                result[i + j] = result[i + j].clone() + ai.clone() * bj.clone();
            }
        }
        DensePoly::new(result)
    }
}
