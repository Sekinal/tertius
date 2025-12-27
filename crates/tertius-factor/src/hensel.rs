//! Hensel lifting for polynomial factorization.
//!
//! Lifts a factorization modulo p to a factorization modulo p^k.
//! Uses quadratic Hensel lifting which doubles precision each iteration.

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

/// Result of Hensel lifting.
#[derive(Clone, Debug)]
pub struct HenselLiftResult {
    /// The lifted factors.
    pub factors: Vec<DensePoly<Z>>,
    /// The prime used.
    pub prime: Integer,
    /// The final modulus (p^k).
    pub modulus: Integer,
    /// Number of lifting steps performed.
    pub steps: usize,
}

/// Lifts a factorization from Z_p to Z_{p^k}.
pub fn hensel_lift(
    f: &DensePoly<Z>,
    factors_mod_p: &[DensePoly<Z>],
    p: &Integer,
    target_k: u32,
) -> HenselLiftResult {
    if factors_mod_p.is_empty() {
        return HenselLiftResult {
            factors: vec![],
            prime: p.clone(),
            modulus: p.clone(),
            steps: 0,
        };
    }

    if factors_mod_p.len() == 1 {
        return HenselLiftResult {
            factors: vec![f.clone()],
            prime: p.clone(),
            modulus: pow_int(p, target_k),
            steps: 0,
        };
    }

    let mut current_factors = factors_mod_p.to_vec();
    let mut modulus = p.clone();
    let mut steps = 0;

    let mut current_k = 1u32;

    while current_k < target_k {
        let next_k = (current_k * 2).min(target_k);
        let new_modulus = pow_int(p, next_k);

        current_factors = lift_factors_one_step(f, &current_factors, &modulus, &new_modulus);

        modulus = new_modulus;
        current_k = next_k;
        steps += 1;
    }

    let half_modulus = modulus.clone() / Integer::new(2);
    let factors: Vec<_> = current_factors
        .into_iter()
        .map(|factor| to_symmetric_rep(&factor, &modulus, &half_modulus))
        .collect();

    HenselLiftResult {
        factors,
        prime: p.clone(),
        modulus,
        steps,
    }
}

fn lift_factors_one_step(
    f: &DensePoly<Z>,
    factors: &[DensePoly<Z>],
    old_modulus: &Integer,
    new_modulus: &Integer,
) -> Vec<DensePoly<Z>> {
    if factors.len() == 2 {
        let (g, h) = (&factors[0], &factors[1]);
        let (new_g, new_h) = hensel_lift_pair(f, g, h, old_modulus, new_modulus);
        vec![new_g, new_h]
    } else {
        let mid = factors.len() / 2;
        let g_factors = &factors[..mid];
        let h_factors = &factors[mid..];

        let g = factors_product(g_factors, new_modulus);
        let h = factors_product(h_factors, new_modulus);

        let (new_g, new_h) = hensel_lift_pair(f, &g, &h, old_modulus, new_modulus);

        let mut result = Vec::with_capacity(factors.len());

        if g_factors.len() == 1 {
            result.push(new_g);
        } else {
            let sub_lifted = lift_factors_one_step(&new_g, g_factors, old_modulus, new_modulus);
            result.extend(sub_lifted);
        }

        if h_factors.len() == 1 {
            result.push(new_h);
        } else {
            let sub_lifted = lift_factors_one_step(&new_h, h_factors, old_modulus, new_modulus);
            result.extend(sub_lifted);
        }

        result
    }
}

fn hensel_lift_pair(
    f: &DensePoly<Z>,
    g: &DensePoly<Z>,
    h: &DensePoly<Z>,
    old_modulus: &Integer,
    new_modulus: &Integer,
) -> (DensePoly<Z>, DensePoly<Z>) {
    let (s, t) = extended_gcd_poly(g, h, old_modulus);

    let gh = poly_mul_mod(g, h, new_modulus);
    let e = poly_sub_mod(f, &gh, new_modulus);

    let te = poly_mul_mod(&t, &e, new_modulus);
    let se = poly_mul_mod(&s, &e, new_modulus);

    let (_, te_mod_g) = poly_div_mod(&te, g, new_modulus);
    let (_, se_mod_h) = poly_div_mod(&se, h, new_modulus);

    let new_g = poly_add_mod(g, &te_mod_g, new_modulus);
    let new_h = poly_add_mod(h, &se_mod_h, new_modulus);

    (new_g, new_h)
}

fn extended_gcd_poly(
    g: &DensePoly<Z>,
    h: &DensePoly<Z>,
    m: &Integer,
) -> (DensePoly<Z>, DensePoly<Z>) {
    let mut old_r = g.clone();
    let mut r = h.clone();
    let mut old_s = DensePoly::new(vec![z_one()]);
    let mut s = DensePoly::zero();
    let mut old_t = DensePoly::zero();
    let mut t = DensePoly::new(vec![z_one()]);

    while !is_zero_mod(&r, m) {
        let (q, new_r) = poly_div_mod(&old_r, &r, m);

        let new_s = poly_sub_mod(&old_s, &poly_mul_mod(&q, &s, m), m);
        let new_t = poly_sub_mod(&old_t, &poly_mul_mod(&q, &t, m), m);

        old_r = r;
        r = new_r;
        old_s = s;
        s = new_s;
        old_t = t;
        t = new_t;
    }

    if !old_r.is_zero() {
        let lead = old_r.coeffs().last().unwrap().clone();
        let lead_inv = mod_inv(&lead.0, m);
        let scale = Z(lead_inv);

        old_s = poly_scale(&old_s, &scale, m);
        old_t = poly_scale(&old_t, &scale, m);
    }

    (old_s, old_t)
}

fn is_zero_mod(p: &DensePoly<Z>, m: &Integer) -> bool {
    p.coeffs().iter().all(|c| (c.0.clone() % m.clone()).is_zero())
}

fn poly_add_mod(a: &DensePoly<Z>, b: &DensePoly<Z>, m: &Integer) -> DensePoly<Z> {
    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let max_len = a_coeffs.len().max(b_coeffs.len());
    let mut result = vec![z_zero(); max_len];

    for (i, c) in a_coeffs.iter().enumerate() {
        result[i] = c.clone();
    }
    for (i, c) in b_coeffs.iter().enumerate() {
        result[i] = Z((result[i].0.clone() + c.0.clone()) % m.clone());
    }

    while result.len() > 1 && result.last().map_or(false, |c| c.0.is_zero()) {
        result.pop();
    }

    DensePoly::new(result)
}

fn poly_sub_mod(a: &DensePoly<Z>, b: &DensePoly<Z>, m: &Integer) -> DensePoly<Z> {
    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let max_len = a_coeffs.len().max(b_coeffs.len());
    let mut result = vec![z_zero(); max_len];

    for (i, c) in a_coeffs.iter().enumerate() {
        result[i] = c.clone();
    }
    for (i, c) in b_coeffs.iter().enumerate() {
        result[i] = Z(((result[i].0.clone() - c.0.clone()) % m.clone() + m.clone()) % m.clone());
    }

    while result.len() > 1 && result.last().map_or(false, |c| c.0.is_zero()) {
        result.pop();
    }

    DensePoly::new(result)
}

fn poly_mul_mod(a: &DensePoly<Z>, b: &DensePoly<Z>, m: &Integer) -> DensePoly<Z> {
    if a.is_zero() || b.is_zero() {
        return DensePoly::zero();
    }

    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let mut result = vec![z_zero(); a_coeffs.len() + b_coeffs.len() - 1];

    for (i, ai) in a_coeffs.iter().enumerate() {
        for (j, bj) in b_coeffs.iter().enumerate() {
            result[i + j] = Z((result[i + j].0.clone() + ai.0.clone() * bj.0.clone()) % m.clone());
        }
    }

    while result.len() > 1 && result.last().map_or(false, |c| c.0.is_zero()) {
        result.pop();
    }

    DensePoly::new(result)
}

fn poly_div_mod(
    a: &DensePoly<Z>,
    b: &DensePoly<Z>,
    m: &Integer,
) -> (DensePoly<Z>, DensePoly<Z>) {
    if b.is_zero() {
        panic!("Division by zero");
    }

    if a.degree() < b.degree() {
        return (DensePoly::zero(), a.clone());
    }

    let mut remainder = a.coeffs().to_vec();
    let divisor = b.coeffs();
    let lead_inv = mod_inv(&divisor.last().unwrap().0, m);
    let deg_diff = a.degree() - b.degree();
    let mut quotient = vec![z_zero(); deg_diff + 1];

    for i in (0..=deg_diff).rev() {
        let idx = i + divisor.len() - 1;
        if idx < remainder.len() {
            let coeff = (remainder[idx].0.clone() * lead_inv.clone()) % m.clone();
            quotient[i] = Z(coeff.clone());

            for (j, d) in divisor.iter().enumerate() {
                let ridx = i + j;
                if ridx < remainder.len() {
                    remainder[ridx] = Z(
                        ((remainder[ridx].0.clone() - coeff.clone() * d.0.clone()) % m.clone()
                            + m.clone())
                            % m.clone(),
                    );
                }
            }
        }
    }

    while remainder.len() > 1 && remainder.last().map_or(false, |c| c.0.is_zero()) {
        remainder.pop();
    }
    while quotient.len() > 1 && quotient.last().map_or(false, |c| c.0.is_zero()) {
        quotient.pop();
    }

    (DensePoly::new(quotient), DensePoly::new(remainder))
}

fn poly_scale(p: &DensePoly<Z>, c: &Z, m: &Integer) -> DensePoly<Z> {
    let coeffs: Vec<_> = p
        .coeffs()
        .iter()
        .map(|coeff| Z((coeff.0.clone() * c.0.clone()) % m.clone()))
        .collect();
    DensePoly::new(coeffs)
}

fn mod_inv(a: &Integer, m: &Integer) -> Integer {
    let mut old_r = a.clone();
    let mut r = m.clone();
    let mut old_s = Integer::new(1);
    let mut s = Integer::new(0);

    while !r.is_zero() {
        let q = old_r.clone() / r.clone();
        let new_r = old_r.clone() - q.clone() * r.clone();
        old_r = r;
        r = new_r;

        let new_s = old_s.clone() - q * s.clone();
        old_s = s;
        s = new_s;
    }

    ((old_s % m.clone()) + m.clone()) % m.clone()
}

fn factors_product(factors: &[DensePoly<Z>], m: &Integer) -> DensePoly<Z> {
    factors
        .iter()
        .fold(DensePoly::new(vec![z_one()]), |acc, f| {
            poly_mul_mod(&acc, f, m)
        })
}

fn to_symmetric_rep(p: &DensePoly<Z>, modulus: &Integer, half: &Integer) -> DensePoly<Z> {
    let coeffs: Vec<_> = p
        .coeffs()
        .iter()
        .map(|c| {
            let c_mod = ((c.0.clone() % modulus.clone()) + modulus.clone()) % modulus.clone();
            if c_mod > *half {
                Z(c_mod - modulus.clone())
            } else {
                Z(c_mod)
            }
        })
        .collect();
    DensePoly::new(coeffs)
}

fn pow_int(p: &Integer, k: u32) -> Integer {
    let mut result = Integer::new(1);
    let mut base = p.clone();
    let mut exp = k;

    while exp > 0 {
        if exp & 1 == 1 {
            result = result * base.clone();
        }
        base = base.clone() * base;
        exp >>= 1;
    }

    result
}

pub fn hensel_lift_batch(
    problems: &[(DensePoly<Z>, Vec<DensePoly<Z>>, Integer, u32)],
) -> Vec<HenselLiftResult> {
    problems
        .par_iter()
        .map(|(f, factors, p, k)| hensel_lift(f, factors, p, *k))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn z(n: i64) -> Z {
        Z(Integer::new(n))
    }

    #[test]
    fn test_mod_inv() {
        let m = Integer::new(7);
        let a = Integer::new(3);
        let inv = mod_inv(&a, &m);
        assert_eq!((a * inv) % m, Integer::new(1));
    }

    #[test]
    fn test_poly_mul_mod() {
        let m = Integer::new(7);
        let a = DensePoly::new(vec![z(1), z(1)]);
        let b = DensePoly::new(vec![z(2), z(1)]);
        let c = poly_mul_mod(&a, &b, &m);

        assert_eq!(c.coeffs()[0].0, Integer::new(2));
        assert_eq!(c.coeffs()[1].0, Integer::new(3));
        assert_eq!(c.coeffs()[2].0, Integer::new(1));
    }

    #[test]
    fn test_hensel_lift_simple() {
        let f = DensePoly::new(vec![z(-1), z(0), z(1)]);
        let g = DensePoly::new(vec![z(4), z(1)]);
        let h = DensePoly::new(vec![z(1), z(1)]);

        let p = Integer::new(5);
        let result = hensel_lift(&f, &[g, h], &p, 2);

        assert_eq!(result.factors.len(), 2);
        assert_eq!(result.modulus, Integer::new(25));
    }

    #[test]
    fn test_pow_int() {
        assert_eq!(pow_int(&Integer::new(2), 0), Integer::new(1));
        assert_eq!(pow_int(&Integer::new(2), 1), Integer::new(2));
        assert_eq!(pow_int(&Integer::new(2), 10), Integer::new(1024));
        assert_eq!(pow_int(&Integer::new(3), 4), Integer::new(81));
    }
}
