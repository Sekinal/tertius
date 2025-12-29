//! Polynomial generation for SIQS.
//!
//! SIQS uses multiple polynomials Q(x) = ax² + bx + c to sieve.
//! We need Q(x) ≡ (ax + b)² (mod n) for the algorithm to work.

use tertius_integers::Integer;

/// A SIQS polynomial Q(x) = ax² + 2bx + c.
#[derive(Debug, Clone)]
pub struct SiqsPolynomial {
    /// The a coefficient (square of product of primes).
    pub a: Integer,
    /// The b coefficient.
    pub b: Integer,
    /// The c coefficient.
    pub c: Integer,
    /// sqrt(a) for computing roots.
    pub sqrt_a: Integer,
    /// The number we're factoring.
    pub n: Integer,
}

impl SiqsPolynomial {
    /// Create a new SIQS polynomial.
    ///
    /// We choose a = q₁²q₂²...qₖ² where qᵢ are primes around √(2n)/M.
    /// Then b is chosen so that b² ≡ n (mod a).
    /// Finally c = (b² - n) / a.
    pub fn new(a: Integer, b: Integer, n: Integer) -> Self {
        let c = (b.clone() * b.clone() - n.clone()) / a.clone();
        let sqrt_a = integer_sqrt(&a);
        Self { a, b, c, sqrt_a, n }
    }

    /// Evaluate Q(x) = ax² + 2bx + c.
    pub fn eval(&self, x: i64) -> Integer {
        let x_int = Integer::from(x);
        let x_sq = x_int.clone() * x_int.clone();
        self.a.clone() * x_sq + Integer::from(2i64) * self.b.clone() * x_int + self.c.clone()
    }

    /// Evaluate (ax + b)² mod n = a * Q(x) + n.
    /// This is what we're actually sieving for.
    pub fn eval_smooth_target(&self, x: i64) -> Integer {
        self.eval(x).abs()
    }
}

/// Integer square root.
fn integer_sqrt(n: &Integer) -> Integer {
    use num_traits::{One, Zero};

    if n.is_zero() {
        return Integer::zero();
    }
    let mut x = n.clone();
    let mut y = (n.clone() + Integer::one()) / Integer::from(2i64);
    while y < x {
        x = y.clone();
        y = (x.clone() + n.clone() / x.clone()) / Integer::from(2i64);
    }
    x
}
