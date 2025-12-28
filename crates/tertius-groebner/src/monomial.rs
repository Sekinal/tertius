//! Packed monomial representation for efficient Gröbner basis computation.
//!
//! Monomials are represented as packed vectors of exponents with efficient
//! comparison, multiplication, and divisibility testing.

use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};

/// Maximum number of variables supported with packed representation.
pub const MAX_VARS: usize = 16;

/// A packed monomial with up to MAX_VARS variables.
///
/// Each exponent is stored as a u16, allowing exponents up to 65535.
/// The total degree is cached for efficient ordering comparisons.
#[derive(Clone, Copy)]
pub struct PackedMonomial {
    /// Exponents for each variable (x_0, x_1, ..., x_{n-1}).
    exponents: [u16; MAX_VARS],
    /// Number of active variables.
    num_vars: u8,
    /// Cached total degree.
    total_degree: u32,
}

impl PackedMonomial {
    /// Creates a new monomial with the given exponents.
    #[must_use]
    pub fn new(exps: &[u16]) -> Self {
        let mut exponents = [0u16; MAX_VARS];
        let n = exps.len().min(MAX_VARS);
        exponents[..n].copy_from_slice(&exps[..n]);

        let total_degree: u32 = exponents.iter().map(|&e| e as u32).sum();

        Self {
            exponents,
            num_vars: n as u8,
            total_degree,
        }
    }

    /// Creates the identity monomial (1).
    #[must_use]
    pub fn one(num_vars: usize) -> Self {
        Self {
            exponents: [0u16; MAX_VARS],
            num_vars: num_vars.min(MAX_VARS) as u8,
            total_degree: 0,
        }
    }

    /// Creates a monomial for a single variable: x_i.
    #[must_use]
    pub fn var(i: usize, num_vars: usize) -> Self {
        let mut exponents = [0u16; MAX_VARS];
        if i < MAX_VARS {
            exponents[i] = 1;
        }
        Self {
            exponents,
            num_vars: num_vars.min(MAX_VARS) as u8,
            total_degree: 1,
        }
    }

    /// Returns the exponent of variable i.
    #[must_use]
    pub fn exponent(&self, i: usize) -> u16 {
        if i < MAX_VARS {
            self.exponents[i]
        } else {
            0
        }
    }

    /// Returns the exponents as a slice.
    #[must_use]
    pub fn exponents(&self) -> &[u16] {
        &self.exponents[..self.num_vars as usize]
    }

    /// Returns the number of variables.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.num_vars as usize
    }

    /// Returns the total degree.
    #[must_use]
    pub fn total_degree(&self) -> u32 {
        self.total_degree
    }

    /// Checks if this is the identity monomial (1).
    #[must_use]
    pub fn is_one(&self) -> bool {
        self.total_degree == 0
    }

    /// Multiplies two monomials.
    #[must_use]
    pub fn mul(&self, other: &Self) -> Self {
        let mut exponents = [0u16; MAX_VARS];
        let n = self.num_vars.max(other.num_vars) as usize;

        for i in 0..n {
            exponents[i] = self.exponents[i].saturating_add(other.exponents[i]);
        }

        Self {
            exponents,
            num_vars: n as u8,
            total_degree: self.total_degree + other.total_degree,
        }
    }

    /// Divides this monomial by another, if divisible.
    ///
    /// Returns `None` if `other` does not divide `self`.
    #[must_use]
    pub fn div(&self, other: &Self) -> Option<Self> {
        if !self.is_divisible_by(other) {
            return None;
        }

        let mut exponents = [0u16; MAX_VARS];
        let n = self.num_vars.max(other.num_vars) as usize;

        for i in 0..n {
            exponents[i] = self.exponents[i] - other.exponents[i];
        }

        Some(Self {
            exponents,
            num_vars: n as u8,
            total_degree: self.total_degree - other.total_degree,
        })
    }

    /// Checks if `other` divides `self`.
    #[must_use]
    pub fn is_divisible_by(&self, other: &Self) -> bool {
        if other.total_degree > self.total_degree {
            return false;
        }

        let n = self.num_vars.max(other.num_vars) as usize;
        for i in 0..n {
            if other.exponents[i] > self.exponents[i] {
                return false;
            }
        }
        true
    }

    /// Computes the least common multiple of two monomials.
    #[must_use]
    pub fn lcm(&self, other: &Self) -> Self {
        let mut exponents = [0u16; MAX_VARS];
        let n = self.num_vars.max(other.num_vars) as usize;

        let mut total = 0u32;
        for i in 0..n {
            exponents[i] = self.exponents[i].max(other.exponents[i]);
            total += exponents[i] as u32;
        }

        Self {
            exponents,
            num_vars: n as u8,
            total_degree: total,
        }
    }

    /// Computes the greatest common divisor of two monomials.
    #[must_use]
    pub fn gcd(&self, other: &Self) -> Self {
        let mut exponents = [0u16; MAX_VARS];
        let n = self.num_vars.max(other.num_vars) as usize;

        let mut total = 0u32;
        for i in 0..n {
            exponents[i] = self.exponents[i].min(other.exponents[i]);
            total += exponents[i] as u32;
        }

        Self {
            exponents,
            num_vars: n as u8,
            total_degree: total,
        }
    }

    /// Checks if two monomials are coprime (GCD = 1).
    #[must_use]
    pub fn is_coprime(&self, other: &Self) -> bool {
        let n = self.num_vars.max(other.num_vars) as usize;
        for i in 0..n {
            if self.exponents[i] > 0 && other.exponents[i] > 0 {
                return false;
            }
        }
        true
    }

    /// Compares using graded reverse lexicographic (grevlex) ordering.
    #[must_use]
    pub fn cmp_grevlex(&self, other: &Self) -> Ordering {
        // First compare total degree
        match self.total_degree.cmp(&other.total_degree) {
            Ordering::Equal => {}
            ord => return ord,
        }

        // Then compare reverse lexicographically (last variable first, reversed)
        let n = self.num_vars.max(other.num_vars) as usize;
        for i in (0..n).rev() {
            match other.exponents[i].cmp(&self.exponents[i]) {
                Ordering::Equal => continue,
                ord => return ord,
            }
        }

        Ordering::Equal
    }

    /// Compares using pure lexicographic ordering.
    #[must_use]
    pub fn cmp_lex(&self, other: &Self) -> Ordering {
        let n = self.num_vars.max(other.num_vars) as usize;
        for i in 0..n {
            match self.exponents[i].cmp(&other.exponents[i]) {
                Ordering::Equal => continue,
                ord => return ord,
            }
        }
        Ordering::Equal
    }
}

impl PartialEq for PackedMonomial {
    fn eq(&self, other: &Self) -> bool {
        if self.total_degree != other.total_degree {
            return false;
        }
        let n = self.num_vars.max(other.num_vars) as usize;
        self.exponents[..n] == other.exponents[..n]
    }
}

impl Eq for PackedMonomial {}

impl Hash for PackedMonomial {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Hash must be consistent with PartialEq which compares up to max(num_vars).
        // To ensure equal monomials hash equally regardless of num_vars,
        // we hash all exponents up to the last non-zero one.
        self.total_degree.hash(state);

        // Find the actual extent of non-zero exponents
        let mut last_nonzero = 0;
        for i in 0..MAX_VARS {
            if self.exponents[i] != 0 {
                last_nonzero = i + 1;
            }
        }
        self.exponents[..last_nonzero].hash(state);
    }
}

impl fmt::Debug for PackedMonomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Mono(")?;
        let mut first = true;
        for (i, &e) in self.exponents[..self.num_vars as usize].iter().enumerate() {
            if e > 0 {
                if !first {
                    write!(f, "*")?;
                }
                first = false;
                if e == 1 {
                    write!(f, "x{}", i)?;
                } else {
                    write!(f, "x{}^{}", i, e)?;
                }
            }
        }
        if first {
            write!(f, "1")?;
        }
        write!(f, ")")
    }
}

impl Default for PackedMonomial {
    fn default() -> Self {
        Self::one(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_monomial_mul() {
        let m1 = PackedMonomial::new(&[1, 2, 0]); // x*y^2
        let m2 = PackedMonomial::new(&[2, 0, 1]); // x^2*z

        let product = m1.mul(&m2);
        assert_eq!(product.exponent(0), 3);
        assert_eq!(product.exponent(1), 2);
        assert_eq!(product.exponent(2), 1);
        assert_eq!(product.total_degree(), 6);
    }

    #[test]
    fn test_monomial_div() {
        let m1 = PackedMonomial::new(&[3, 2, 1]); // x^3*y^2*z
        let m2 = PackedMonomial::new(&[1, 1, 0]); // x*y

        let quotient = m1.div(&m2).unwrap();
        assert_eq!(quotient.exponent(0), 2);
        assert_eq!(quotient.exponent(1), 1);
        assert_eq!(quotient.exponent(2), 1);

        // Can't divide x by x^2
        let m3 = PackedMonomial::new(&[1, 0, 0]);
        let m4 = PackedMonomial::new(&[2, 0, 0]);
        assert!(m3.div(&m4).is_none());
    }

    #[test]
    fn test_monomial_lcm_gcd() {
        let m1 = PackedMonomial::new(&[2, 1, 0]); // x^2*y
        let m2 = PackedMonomial::new(&[1, 3, 0]); // x*y^3

        let lcm = m1.lcm(&m2);
        assert_eq!(lcm.exponent(0), 2);
        assert_eq!(lcm.exponent(1), 3);

        let gcd = m1.gcd(&m2);
        assert_eq!(gcd.exponent(0), 1);
        assert_eq!(gcd.exponent(1), 1);
    }

    #[test]
    fn test_grevlex_ordering() {
        // x^2 > xy > y^2 > x > y > 1 in grevlex
        let x2 = PackedMonomial::new(&[2, 0]);
        let xy = PackedMonomial::new(&[1, 1]);
        let y2 = PackedMonomial::new(&[0, 2]);
        let x = PackedMonomial::new(&[1, 0]);
        let y = PackedMonomial::new(&[0, 1]);
        let one = PackedMonomial::one(2);

        assert_eq!(x2.cmp_grevlex(&xy), Ordering::Greater);
        assert_eq!(xy.cmp_grevlex(&y2), Ordering::Greater);
        assert_eq!(y2.cmp_grevlex(&x), Ordering::Greater);
        assert_eq!(x.cmp_grevlex(&y), Ordering::Greater);
        assert_eq!(y.cmp_grevlex(&one), Ordering::Greater);
    }

    #[test]
    fn test_coprime() {
        let m1 = PackedMonomial::new(&[2, 0, 0]); // x^2
        let m2 = PackedMonomial::new(&[0, 1, 0]); // y

        assert!(m1.is_coprime(&m2));

        let m3 = PackedMonomial::new(&[1, 1, 0]); // xy
        assert!(!m1.is_coprime(&m3));
    }
}
