//! F5 signatures for Gröbner basis computation.
//!
//! Signatures track the "origin" of polynomials through reductions,
//! enabling early detection of useless S-polynomials.

use crate::monomial::PackedMonomial;
use std::cmp::Ordering;

/// An F5 signature represents the "origin" of a polynomial.
///
/// A signature (m, i) means this polynomial was derived from multiplying
/// generator g_i by monomial m, after reductions.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct Signature {
    /// The monomial multiplier.
    pub monomial: PackedMonomial,
    /// Index of the original generator (0-indexed).
    pub index: usize,
}

impl Signature {
    /// Creates a new signature.
    #[must_use]
    pub fn new(monomial: PackedMonomial, index: usize) -> Self {
        Self { monomial, index }
    }

    /// Creates the signature for a generator g_i.
    #[must_use]
    pub fn generator(index: usize, num_vars: usize) -> Self {
        Self {
            monomial: PackedMonomial::one(num_vars),
            index,
        }
    }

    /// Multiplies this signature by a monomial.
    #[must_use]
    pub fn mul(&self, m: &PackedMonomial) -> Self {
        Self {
            monomial: self.monomial.mul(m),
            index: self.index,
        }
    }

    /// Compares signatures using Position-Over-Term (POT) ordering.
    ///
    /// POT ordering: (m1, i1) < (m2, i2) iff i1 < i2, or i1 == i2 and m1 < m2.
    #[must_use]
    pub fn cmp_pot(&self, other: &Self) -> Ordering {
        match self.index.cmp(&other.index) {
            Ordering::Equal => self.monomial.cmp_grevlex(&other.monomial),
            ord => ord,
        }
    }

    /// Compares signatures using Term-Over-Position (TOP) ordering.
    ///
    /// TOP ordering: (m1, i1) < (m2, i2) iff m1 < m2, or m1 == m2 and i1 < i2.
    #[must_use]
    pub fn cmp_top(&self, other: &Self) -> Ordering {
        match self.monomial.cmp_grevlex(&other.monomial) {
            Ordering::Equal => self.index.cmp(&other.index),
            ord => ord,
        }
    }

    /// Checks if this signature divides another.
    #[must_use]
    pub fn divides(&self, other: &Self) -> bool {
        self.index == other.index && other.monomial.is_divisible_by(&self.monomial)
    }
}

impl Default for Signature {
    fn default() -> Self {
        Self::generator(0, 0)
    }
}

impl Ord for Signature {
    fn cmp(&self, other: &Self) -> Ordering {
        // Default to POT ordering
        self.cmp_pot(other)
    }
}

impl PartialOrd for Signature {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// A critical pair with signatures for the S-polynomial.
#[derive(Clone, Debug)]
pub struct SignedPair {
    /// Index of first polynomial in basis.
    pub i: usize,
    /// Index of second polynomial in basis.
    pub j: usize,
    /// Signature of the S-polynomial.
    pub signature: Signature,
    /// The LCM of the leading monomials.
    pub lcm: PackedMonomial,
    /// Sugar degree (for selection strategy).
    pub sugar: u32,
}

impl SignedPair {
    /// Creates a new signed pair.
    #[must_use]
    pub fn new(
        i: usize,
        j: usize,
        sig_i: &Signature,
        sig_j: &Signature,
        lm_i: &PackedMonomial,
        lm_j: &PackedMonomial,
        sugar_i: u32,
        sugar_j: u32,
    ) -> Self {
        let lcm = lm_i.lcm(lm_j);

        // Compute the multipliers for each polynomial
        let mult_i = lcm.div(lm_i).unwrap();
        let mult_j = lcm.div(lm_j).unwrap();

        // The signature of S(f_i, f_j) is max(mult_i * sig_i, mult_j * sig_j)
        let sig_spoly_i = sig_i.mul(&mult_i);
        let sig_spoly_j = sig_j.mul(&mult_j);

        let signature = if sig_spoly_i.cmp_pot(&sig_spoly_j) == Ordering::Greater {
            sig_spoly_i
        } else {
            sig_spoly_j
        };

        // Sugar degree is max of (sugar_i + deg(mult_i), sugar_j + deg(mult_j))
        let sugar = (sugar_i + mult_i.total_degree()).max(sugar_j + mult_j.total_degree());

        Self {
            i,
            j,
            signature,
            lcm,
            sugar,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_signature_pot_ordering() {
        let m1 = PackedMonomial::new(&[1, 0]);
        let m2 = PackedMonomial::new(&[0, 1]);

        let sig1 = Signature::new(m1, 0);
        let sig2 = Signature::new(m2, 1);

        // sig1 < sig2 because index 0 < index 1
        assert_eq!(sig1.cmp_pot(&sig2), Ordering::Less);

        let sig3 = Signature::new(m2, 0);
        // sig1 > sig3 because same index, and x > y in grevlex
        assert_eq!(sig1.cmp_pot(&sig3), Ordering::Greater);
    }

    #[test]
    fn test_signature_mul() {
        let m = PackedMonomial::new(&[1, 1]);
        let sig = Signature::generator(0, 2);

        let product = sig.mul(&m);
        assert_eq!(product.index, 0);
        assert_eq!(product.monomial.exponent(0), 1);
        assert_eq!(product.monomial.exponent(1), 1);
    }

    #[test]
    fn test_signature_divides() {
        let sig1 = Signature::new(PackedMonomial::new(&[1, 0]), 0);
        let sig2 = Signature::new(PackedMonomial::new(&[2, 1]), 0);
        let sig3 = Signature::new(PackedMonomial::new(&[1, 0]), 1);

        assert!(sig1.divides(&sig2)); // x divides x^2*y
        assert!(!sig2.divides(&sig1)); // x^2*y doesn't divide x
        assert!(!sig1.divides(&sig3)); // Different index
    }
}
