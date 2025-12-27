//! Squarefree factorization of polynomials.
//!
//! Produces a factorization where each factor is squarefree
//! (no repeated roots).

use num_traits::{One, Zero};
use rayon::prelude::*;
use tertius_integers::Integer;
use tertius_poly::dense::DensePoly;
use tertius_rings::integers::Z;
use tertius_rings::traits::Ring;

// Helper functions to avoid trait ambiguity
fn z_zero() -> Z {
    <Z as Ring>::zero()
}

fn z_one() -> Z {
    <Z as Ring>::one()
}

/// A squarefree factor with its multiplicity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SquarefreeFactor {
    /// The squarefree polynomial.
    pub factor: DensePoly<Z>,
    /// The multiplicity (power) of this factor.
    pub multiplicity: usize,
}

/// Result of squarefree factorization.
#[derive(Clone, Debug)]
pub struct SquarefreeFactorization {
    /// Leading coefficient (content that was factored out).
    pub content: Z,
    /// List of squarefree factors with multiplicities.
    pub factors: Vec<SquarefreeFactor>,
}

impl SquarefreeFactorization {
    /// Reconstructs the original polynomial from factorization.
    pub fn to_polynomial(&self) -> DensePoly<Z> {
        let mut result = DensePoly::new(vec![self.content.clone()]);

        for sf in &self.factors {
            for _ in 0..sf.multiplicity {
                result = poly_mul(&result, &sf.factor);
            }
        }

        result
    }
}

/// Computes the squarefree factorization of a polynomial over Z.
pub fn squarefree_factorization(f: &DensePoly<Z>) -> SquarefreeFactorization {
    if f.is_zero() {
        return SquarefreeFactorization {
            content: z_zero(),
            factors: vec![],
        };
    }

    if f.degree() == 0 {
        return SquarefreeFactorization {
            content: f.coeffs()[0].clone(),
            factors: vec![],
        };
    }

    let content = polynomial_content(f);
    let primitive = if content.0.is_one() {
        f.clone()
    } else {
        poly_div_exact(f, &DensePoly::new(vec![content.clone()]))
    };

    let (sign, primitive) = if primitive.coeffs().last().unwrap().0 < Integer::new(0) {
        (Z(Integer::new(-1)), poly_negate(&primitive))
    } else {
        (z_one(), primitive)
    };

    let final_content = Z(content.0 * sign.0);

    let factors = yun_algorithm(&primitive);

    SquarefreeFactorization {
        content: final_content,
        factors,
    }
}

fn yun_algorithm(f: &DensePoly<Z>) -> Vec<SquarefreeFactor> {
    let mut factors = Vec::new();

    let f_prime = polynomial_derivative(f);

    if f_prime.is_zero() {
        if f.degree() > 0 {
            factors.push(SquarefreeFactor {
                factor: f.clone(),
                multiplicity: 1,
            });
        }
        return factors;
    }

    let a0 = poly_gcd(f, &f_prime);
    let b1 = poly_div_exact(f, &a0);
    let c1 = poly_div_exact(&f_prime, &a0);
    let b1_prime = polynomial_derivative(&b1);
    let d1 = poly_sub(&c1, &b1_prime);

    let mut b = b1;
    let mut d = d1;
    let mut i = 1;

    while b.degree() > 0 {
        let a = poly_gcd(&b, &d);

        if a.degree() > 0 {
            factors.push(SquarefreeFactor {
                factor: make_primitive(&a),
                multiplicity: i,
            });
        }

        let new_b = poly_div_exact(&b, &a);
        let c = poly_div_exact(&d, &a);
        let new_b_prime = polynomial_derivative(&new_b);
        let new_d = poly_sub(&c, &new_b_prime);

        b = new_b;
        d = new_d;
        i += 1;
    }

    factors
        .into_iter()
        .filter(|sf| sf.factor.degree() > 0)
        .collect()
}

fn polynomial_derivative(f: &DensePoly<Z>) -> DensePoly<Z> {
    let coeffs = f.coeffs();
    if coeffs.len() <= 1 {
        return DensePoly::zero();
    }

    let deriv: Vec<Z> = coeffs
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| Z(c.0.clone() * Integer::new(i as i64)))
        .collect();

    DensePoly::new(deriv)
}

fn polynomial_content(f: &DensePoly<Z>) -> Z {
    f.coeffs()
        .iter()
        .fold(z_zero(), |acc, c| Z(gcd_int(&acc.0, &c.0)))
}

fn gcd_int(a: &Integer, b: &Integer) -> Integer {
    if b.is_zero() {
        if *a < Integer::new(0) {
            -a.clone()
        } else {
            a.clone()
        }
    } else {
        gcd_int(b, &(a.clone() % b.clone()))
    }
}

fn poly_gcd(a: &DensePoly<Z>, b: &DensePoly<Z>) -> DensePoly<Z> {
    if b.is_zero() {
        return make_primitive(a);
    }
    if a.is_zero() {
        return make_primitive(b);
    }

    let (mut a, mut b) = if a.degree() >= b.degree() {
        (a.clone(), b.clone())
    } else {
        (b.clone(), a.clone())
    };

    let mut g = z_one();
    let mut h = z_one();

    while !b.is_zero() && b.degree() > 0 {
        let delta = a.degree() - b.degree();

        let r = pseudo_remainder(&a, &b);

        if r.is_zero() || r.degree() == 0 {
            break;
        }

        let g_h_delta = Z(g.0.clone() * pow_z(&h.0, delta as u32));
        let r_normalized = poly_div_content(&r, &g_h_delta);

        a = b;
        b = r_normalized;
        g = Z(a.coeffs().last().unwrap().0.clone().abs());

        if delta == 0 {
            // h stays the same
        } else if delta == 1 {
            h = g.clone();
        } else {
            h = Z(pow_z(&g.0, delta as u32) / pow_z(&h.0, (delta - 1) as u32));
        }
    }

    make_primitive(&a)
}

fn pseudo_remainder(a: &DensePoly<Z>, b: &DensePoly<Z>) -> DensePoly<Z> {
    if b.is_zero() {
        panic!("Division by zero");
    }

    let deg_a = a.degree();
    let deg_b = b.degree();

    if deg_a < deg_b {
        return a.clone();
    }

    let mut r = a.coeffs().to_vec();
    let b_coeffs = b.coeffs();
    let lead_b = b_coeffs.last().unwrap();
    let delta = deg_a - deg_b;

    let scale = pow_z(&lead_b.0, (delta + 1) as u32);
    for c in &mut r {
        *c = Z(c.0.clone() * scale.clone());
    }

    for i in (0..=delta).rev() {
        let idx = i + deg_b;
        if idx < r.len() {
            let q = r[idx].0.clone() / lead_b.0.clone();
            for (j, bc) in b_coeffs.iter().enumerate() {
                let ridx = i + j;
                r[ridx] = Z(r[ridx].0.clone() - q.clone() * bc.0.clone());
            }
        }
    }

    while r.len() > 1 && r.last().map_or(false, |c| c.0.is_zero()) {
        r.pop();
    }

    DensePoly::new(r)
}

fn pow_z(base: &Integer, exp: u32) -> Integer {
    let mut result = Integer::new(1);
    let mut b = base.clone();
    let mut e = exp;

    while e > 0 {
        if e & 1 == 1 {
            result = result * b.clone();
        }
        b = b.clone() * b;
        e >>= 1;
    }

    result
}

fn poly_div_content(f: &DensePoly<Z>, c: &Z) -> DensePoly<Z> {
    let coeffs: Vec<Z> = f
        .coeffs()
        .iter()
        .map(|coeff| Z(coeff.0.clone() / c.0.clone()))
        .collect();
    DensePoly::new(coeffs)
}

fn make_primitive(f: &DensePoly<Z>) -> DensePoly<Z> {
    if f.is_zero() {
        return f.clone();
    }

    let content = polynomial_content(f);
    if content.0.is_one() {
        let lead = f.coeffs().last().unwrap();
        if lead.0 < Integer::new(0) {
            poly_negate(f)
        } else {
            f.clone()
        }
    } else {
        let mut result = poly_div_content(f, &content);
        if result.coeffs().last().unwrap().0 < Integer::new(0) {
            result = poly_negate(&result);
        }
        result
    }
}

fn poly_div_exact(a: &DensePoly<Z>, b: &DensePoly<Z>) -> DensePoly<Z> {
    if b.is_zero() {
        panic!("Division by zero");
    }

    if a.degree() < b.degree() {
        return DensePoly::zero();
    }

    if b.degree() == 0 {
        let c = &b.coeffs()[0];
        return poly_div_content(a, c);
    }

    let mut remainder = a.coeffs().to_vec();
    let divisor = b.coeffs();
    let deg_diff = a.degree() - b.degree();
    let mut quotient = vec![z_zero(); deg_diff + 1];

    let lead_b = divisor.last().unwrap();

    for i in (0..=deg_diff).rev() {
        let idx = i + divisor.len() - 1;
        if idx < remainder.len() && !remainder[idx].0.is_zero() {
            let q = remainder[idx].0.clone() / lead_b.0.clone();
            let rem = remainder[idx].0.clone() % lead_b.0.clone();
            if !rem.is_zero() {
                continue;
            }

            quotient[i] = Z(q.clone());

            for (j, d) in divisor.iter().enumerate() {
                let ridx = i + j;
                remainder[ridx] = Z(remainder[ridx].0.clone() - q.clone() * d.0.clone());
            }
        }
    }

    while quotient.len() > 1 && quotient.last().map_or(false, |c| c.0.is_zero()) {
        quotient.pop();
    }

    DensePoly::new(quotient)
}

fn poly_mul(a: &DensePoly<Z>, b: &DensePoly<Z>) -> DensePoly<Z> {
    if a.is_zero() || b.is_zero() {
        return DensePoly::zero();
    }

    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let mut result = vec![z_zero(); a_coeffs.len() + b_coeffs.len() - 1];

    for (i, ai) in a_coeffs.iter().enumerate() {
        for (j, bj) in b_coeffs.iter().enumerate() {
            result[i + j] = Z(result[i + j].0.clone() + ai.0.clone() * bj.0.clone());
        }
    }

    DensePoly::new(result)
}

fn poly_sub(a: &DensePoly<Z>, b: &DensePoly<Z>) -> DensePoly<Z> {
    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let max_len = a_coeffs.len().max(b_coeffs.len());
    let mut result = vec![z_zero(); max_len];

    for (i, c) in a_coeffs.iter().enumerate() {
        result[i] = c.clone();
    }
    for (i, c) in b_coeffs.iter().enumerate() {
        result[i] = Z(result[i].0.clone() - c.0.clone());
    }

    while result.len() > 1 && result.last().map_or(false, |c| c.0.is_zero()) {
        result.pop();
    }

    DensePoly::new(result)
}

fn poly_negate(f: &DensePoly<Z>) -> DensePoly<Z> {
    let coeffs: Vec<Z> = f.coeffs().iter().map(|c| Z(-c.0.clone())).collect();
    DensePoly::new(coeffs)
}

pub fn squarefree_factorization_batch(polys: &[DensePoly<Z>]) -> Vec<SquarefreeFactorization> {
    polys.par_iter().map(squarefree_factorization).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn z(n: i64) -> Z {
        Z(Integer::new(n))
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Z> {
        DensePoly::new(coeffs.iter().map(|&n| z(n)).collect())
    }

    #[test]
    fn test_derivative() {
        let f = poly(&[1, 1, 2, 1]);
        let df = polynomial_derivative(&f);

        assert_eq!(df.coeffs().len(), 3);
        assert_eq!(df.coeffs()[0].0, Integer::new(1));
        assert_eq!(df.coeffs()[1].0, Integer::new(4));
        assert_eq!(df.coeffs()[2].0, Integer::new(3));
    }

    #[test]
    fn test_content() {
        let f = poly(&[6, 12, 18]);
        let c = polynomial_content(&f);
        assert_eq!(c.0, Integer::new(6));
    }

    #[test]
    fn test_gcd_int() {
        assert_eq!(gcd_int(&Integer::new(12), &Integer::new(8)), Integer::new(4));
        assert_eq!(gcd_int(&Integer::new(17), &Integer::new(5)), Integer::new(1));
    }

    #[test]
    fn test_squarefree_simple() {
        let f = poly(&[1, 2, 1]);
        let sf = squarefree_factorization(&f);

        assert_eq!(sf.factors.len(), 1);
        assert_eq!(sf.factors[0].multiplicity, 2);
        assert_eq!(sf.factors[0].factor.degree(), 1);
    }

    #[test]
    fn test_squarefree_mixed() {
        let xp1 = poly(&[1, 1]);
        let xp2 = poly(&[2, 1]);
        let xp1_sq = poly_mul(&xp1, &xp1);
        let f = poly_mul(&xp1_sq, &xp2);

        let sf = squarefree_factorization(&f);

        assert_eq!(sf.factors.len(), 2);

        let total_degree: usize = sf.factors.iter().map(|f| f.factor.degree() * f.multiplicity).sum();
        assert_eq!(total_degree, 3);
    }

    #[test]
    fn test_squarefree_already() {
        let f = poly(&[2, 3, 1]);
        let sf = squarefree_factorization(&f);

        assert_eq!(sf.factors.len(), 1);
        assert_eq!(sf.factors[0].multiplicity, 1);
    }

    #[test]
    fn test_reconstruction() {
        let f = poly(&[1, 2, 1]);
        let sf = squarefree_factorization(&f);
        let reconstructed = sf.to_polynomial();

        assert_eq!(f.coeffs(), reconstructed.coeffs());
    }
}
