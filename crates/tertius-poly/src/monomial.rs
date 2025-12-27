//! Bit-packed monomials for efficient multivariate polynomial arithmetic.
//!
//! Monomials are packed into a single u64 or u128 for fast comparison
//! and multiplication using SIMD-friendly integer operations.

use std::cmp::Ordering;

/// A bit-packed monomial representation.
///
/// For up to 3 variables, each exponent gets 21 bits (max degree 2^21 - 1).
/// For more variables, a dynamic representation is used.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub struct PackedMonomial(u64);

impl PackedMonomial {
    /// Bits per exponent for small variable counts.
    const BITS_PER_EXP: u32 = 21;
    /// Maximum variables that fit in u64.
    const MAX_VARS_U64: usize = 3;
    /// Mask for a single exponent.
    const EXP_MASK: u64 = (1 << Self::BITS_PER_EXP) - 1;

    /// Creates the monomial 1 (all exponents zero).
    #[must_use]
    pub const fn one(_num_vars: usize) -> Self {
        Self(0)
    }

    /// Creates the monomial x_i.
    #[must_use]
    pub fn var(i: usize, num_vars: usize) -> Self {
        assert!(i < num_vars);
        assert!(num_vars <= Self::MAX_VARS_U64);
        Self(1 << (i as u32 * Self::BITS_PER_EXP))
    }

    /// Creates a monomial from exponents.
    #[must_use]
    pub fn from_exponents(exps: &[u32]) -> Self {
        assert!(exps.len() <= Self::MAX_VARS_U64);

        let mut packed = 0u64;
        for (i, &e) in exps.iter().enumerate() {
            assert!(e <= Self::EXP_MASK as u32);
            packed |= (e as u64) << (i as u32 * Self::BITS_PER_EXP);
        }

        Self(packed)
    }

    /// Returns the exponent of variable i.
    #[must_use]
    pub fn exponent(&self, i: usize) -> u32 {
        ((self.0 >> (i as u32 * Self::BITS_PER_EXP)) & Self::EXP_MASK) as u32
    }

    /// Returns all exponents.
    #[must_use]
    pub fn exponents(&self, num_vars: usize) -> Vec<u32> {
        (0..num_vars).map(|i| self.exponent(i)).collect()
    }

    /// Multiplies two monomials (adds exponents).
    #[must_use]
    pub fn mul(&self, other: &Self) -> Self {
        // Exponent addition is just integer addition when packed
        Self(self.0 + other.0)
    }

    /// Divides by another monomial if possible.
    ///
    /// Returns `Some(quotient)` if self is divisible by other.
    #[must_use]
    pub fn div(&self, other: &Self, num_vars: usize) -> Option<Self> {
        // Check each exponent
        for i in 0..num_vars {
            if self.exponent(i) < other.exponent(i) {
                return None;
            }
        }
        Some(Self(self.0 - other.0))
    }

    /// Returns true if self is divisible by other.
    #[must_use]
    pub fn divides(&self, other: &Self, num_vars: usize) -> bool {
        for i in 0..num_vars {
            if self.exponent(i) > other.exponent(i) {
                return false;
            }
        }
        true
    }

    /// Computes the total degree.
    #[must_use]
    pub fn total_degree(&self, num_vars: usize) -> u32 {
        (0..num_vars).map(|i| self.exponent(i)).sum()
    }

    /// Computes the least common multiple of two monomials.
    #[must_use]
    pub fn lcm(&self, other: &Self, num_vars: usize) -> Self {
        let mut result = 0u64;
        for i in 0..num_vars {
            let e = self.exponent(i).max(other.exponent(i));
            result |= (e as u64) << (i as u32 * Self::BITS_PER_EXP);
        }
        Self(result)
    }

    /// Computes the greatest common divisor of two monomials.
    #[must_use]
    pub fn gcd(&self, other: &Self, num_vars: usize) -> Self {
        let mut result = 0u64;
        for i in 0..num_vars {
            let e = self.exponent(i).min(other.exponent(i));
            result |= (e as u64) << (i as u32 * Self::BITS_PER_EXP);
        }
        Self(result)
    }

    /// Returns the raw packed value.
    #[must_use]
    pub const fn raw(&self) -> u64 {
        self.0
    }

    /// Converts to a human-readable string.
    #[must_use]
    pub fn to_string(&self, num_vars: usize) -> String {
        let vars = ['x', 'y', 'z', 'w', 'u', 'v'];
        let mut parts = Vec::new();

        for i in 0..num_vars {
            let e = self.exponent(i);
            if e > 0 {
                let var_name = if i < vars.len() {
                    vars[i].to_string()
                } else {
                    format!("x{i}")
                };

                if e == 1 {
                    parts.push(var_name);
                } else {
                    parts.push(format!("{var_name}^{e}"));
                }
            }
        }

        if parts.is_empty() {
            "1".to_string()
        } else {
            parts.join("*")
        }
    }
}

impl std::fmt::Display for PackedMonomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Default to 3 variables
        write!(f, "{}", self.to_string(3))
    }
}

/// Compares two monomials lexicographically.
pub fn cmp_lex(a: &PackedMonomial, b: &PackedMonomial, num_vars: usize) -> Ordering {
    for i in 0..num_vars {
        match a.exponent(i).cmp(&b.exponent(i)) {
            Ordering::Equal => continue,
            ord => return ord,
        }
    }
    Ordering::Equal
}

/// Compares two monomials by graded reverse lexicographic order.
pub fn cmp_grevlex(a: &PackedMonomial, b: &PackedMonomial, num_vars: usize) -> Ordering {
    // First compare total degree
    match a.total_degree(num_vars).cmp(&b.total_degree(num_vars)) {
        Ordering::Equal => {}
        ord => return ord,
    }

    // Then compare in reverse order, reversed
    for i in (0..num_vars).rev() {
        match b.exponent(i).cmp(&a.exponent(i)) {
            Ordering::Equal => continue,
            ord => return ord,
        }
    }
    Ordering::Equal
}

/// Compares two monomials by graded lexicographic order.
pub fn cmp_grlex(a: &PackedMonomial, b: &PackedMonomial, num_vars: usize) -> Ordering {
    // First compare total degree
    match a.total_degree(num_vars).cmp(&b.total_degree(num_vars)) {
        Ordering::Equal => {}
        ord => return ord,
    }

    // Then lexicographic
    cmp_lex(a, b, num_vars)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic() {
        let x = PackedMonomial::var(0, 3);
        let y = PackedMonomial::var(1, 3);

        assert_eq!(x.exponent(0), 1);
        assert_eq!(x.exponent(1), 0);
        assert_eq!(y.exponent(0), 0);
        assert_eq!(y.exponent(1), 1);
    }

    #[test]
    fn test_mul() {
        let x = PackedMonomial::var(0, 3);
        let y = PackedMonomial::var(1, 3);

        let xy = x.mul(&y);
        assert_eq!(xy.exponent(0), 1);
        assert_eq!(xy.exponent(1), 1);

        let x2y = x.mul(&xy);
        assert_eq!(x2y.exponent(0), 2);
        assert_eq!(x2y.exponent(1), 1);
    }

    #[test]
    fn test_div() {
        let x2y = PackedMonomial::from_exponents(&[2, 1, 0]);
        let xy = PackedMonomial::from_exponents(&[1, 1, 0]);
        let x = PackedMonomial::var(0, 3);

        assert_eq!(x2y.div(&xy, 3), Some(x));
        assert_eq!(xy.div(&x2y, 3), None);
    }

    #[test]
    fn test_grevlex_order() {
        let x2 = PackedMonomial::from_exponents(&[2, 0]);
        let xy = PackedMonomial::from_exponents(&[1, 1]);
        let y2 = PackedMonomial::from_exponents(&[0, 2]);

        // In grevlex: x^2 > xy > y^2
        assert_eq!(cmp_grevlex(&x2, &xy, 2), Ordering::Greater);
        assert_eq!(cmp_grevlex(&xy, &y2, 2), Ordering::Greater);
    }
}
