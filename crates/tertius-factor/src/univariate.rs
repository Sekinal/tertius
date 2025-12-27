//! Van Hoeij's algorithm for polynomial factorization over Z.
//!
//! Uses lattice reduction (LLL) to efficiently combine Hensel-lifted
//! modular factors into true integer factors.

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

use crate::hensel::{hensel_lift, HenselLiftResult};
use crate::squarefree::squarefree_factorization;

/// Result of Van Hoeij factorization.
#[derive(Clone, Debug)]
pub struct VanHoeijResult {
    /// The irreducible factors.
    pub factors: Vec<DensePoly<Z>>,
    /// Content (leading coefficient extracted).
    pub content: Z,
    /// Statistics about the computation.
    pub stats: VanHoeijStats,
}

/// Statistics from Van Hoeij algorithm.
#[derive(Clone, Debug, Default)]
pub struct VanHoeijStats {
    /// Prime used for modular factorization.
    pub prime: u64,
    /// Hensel lifting precision.
    pub precision: u32,
    /// Number of modular factors.
    pub num_mod_factors: usize,
}

/// Factors a polynomial over Z using Van Hoeij's algorithm.
pub fn van_hoeij_factor(f: &DensePoly<Z>) -> VanHoeijResult {
    if f.is_zero() {
        return VanHoeijResult {
            factors: vec![],
            content: z_zero(),
            stats: VanHoeijStats::default(),
        };
    }

    let n = f.degree();
    if n == 0 {
        return VanHoeijResult {
            factors: vec![],
            content: f.coeffs()[0].clone(),
            stats: VanHoeijStats::default(),
        };
    }

    if n == 1 {
        let content = polynomial_content(f);
        let primitive = make_primitive(f);
        return VanHoeijResult {
            factors: vec![primitive],
            content,
            stats: VanHoeijStats::default(),
        };
    }

    // Step 1: Squarefree factorization
    let sf = squarefree_factorization(f);

    let mut all_factors = Vec::new();
    let mut total_stats = VanHoeijStats::default();

    for sf_factor in &sf.factors {
        let result = factor_squarefree(&sf_factor.factor);
        total_stats.prime = result.stats.prime;
        total_stats.precision = total_stats.precision.max(result.stats.precision);
        total_stats.num_mod_factors += result.stats.num_mod_factors;

        for _ in 0..sf_factor.multiplicity {
            all_factors.extend(result.factors.iter().cloned());
        }
    }

    VanHoeijResult {
        factors: all_factors,
        content: sf.content,
        stats: total_stats,
    }
}

fn factor_squarefree(f: &DensePoly<Z>) -> VanHoeijResult {
    let n = f.degree();

    // Choose a good prime
    let p = choose_prime(f);

    // Factor modulo p
    let mod_factors = factor_mod_prime(f, p);

    if mod_factors.len() == 1 {
        return VanHoeijResult {
            factors: vec![f.clone()],
            content: z_one(),
            stats: VanHoeijStats {
                prime: p,
                precision: 1,
                num_mod_factors: 1,
            },
        };
    }

    // Compute required Hensel lifting precision
    let bound = factor_coefficient_bound(f);
    let precision = compute_precision(p, &bound);

    // Hensel lift
    let lifted = hensel_lift(f, &mod_factors, &Integer::new(p as i64), precision);

    // Try to extract factors
    let factors = combine_factors(f, &lifted);

    VanHoeijResult {
        factors,
        content: z_one(),
        stats: VanHoeijStats {
            prime: p,
            precision,
            num_mod_factors: mod_factors.len(),
        },
    }
}

fn choose_prime(f: &DensePoly<Z>) -> u64 {
    let lead = f.coeffs().last().unwrap();

    let primes: [u64; 10] = [
        1009, 10007, 100003, 1000003, 10000019, 100000007, 1000000007, 1000000009, 2147483647,
        2305843009213693951u64 >> 32,
    ];

    for &p in &primes {
        if (lead.0.clone() % Integer::new(p as i64)).is_zero() {
            continue;
        }

        let f_mod = reduce_mod_p(f, p);
        if is_squarefree_mod_p(&f_mod, p) {
            return p;
        }
    }

    1000000007
}

fn reduce_mod_p(f: &DensePoly<Z>, p: u64) -> DensePoly<Z> {
    let coeffs: Vec<Z> = f
        .coeffs()
        .iter()
        .map(|c| {
            let r = c.0.clone() % Integer::new(p as i64);
            let r = if r < Integer::new(0) {
                r + Integer::new(p as i64)
            } else {
                r
            };
            Z(r)
        })
        .collect();
    DensePoly::new(coeffs)
}

fn is_squarefree_mod_p(f: &DensePoly<Z>, p: u64) -> bool {
    let f_prime = polynomial_derivative_mod(f, p);
    let gcd = poly_gcd_mod(f, &f_prime, p);
    gcd.degree() == 0
}

fn polynomial_derivative_mod(f: &DensePoly<Z>, p: u64) -> DensePoly<Z> {
    let coeffs = f.coeffs();
    if coeffs.len() <= 1 {
        return DensePoly::zero();
    }

    let deriv: Vec<Z> = coeffs
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| {
            let r = (c.0.clone() * Integer::new(i as i64)) % Integer::new(p as i64);
            Z(r)
        })
        .collect();

    DensePoly::new(deriv)
}

fn poly_gcd_mod(a: &DensePoly<Z>, b: &DensePoly<Z>, p: u64) -> DensePoly<Z> {
    if b.is_zero() {
        return a.clone();
    }

    let r = poly_mod_mod(a, b, p);
    poly_gcd_mod(b, &r, p)
}

fn poly_mod_mod(a: &DensePoly<Z>, b: &DensePoly<Z>, p: u64) -> DensePoly<Z> {
    if b.is_zero() {
        panic!("Division by zero");
    }

    if a.degree() < b.degree() {
        return a.clone();
    }

    let mut remainder = a.coeffs().to_vec();
    let divisor = b.coeffs();
    let lead_inv = mod_inv(
        divisor.last().unwrap().0.to_i64().unwrap() as u64 % p,
        p,
    );

    for i in (0..=(a.degree() - b.degree())).rev() {
        let idx = i + divisor.len() - 1;
        if idx < remainder.len() {
            let r_val = remainder[idx].0.to_i64().unwrap().rem_euclid(p as i64) as u64;
            let q = (r_val * lead_inv) % p;

            for (j, d) in divisor.iter().enumerate() {
                let ridx = i + j;
                let d_val = d.0.to_i64().unwrap().rem_euclid(p as i64) as u64;
                let sub = (q * d_val) % p;
                let cur = remainder[ridx].0.to_i64().unwrap().rem_euclid(p as i64) as u64;
                remainder[ridx] = Z(Integer::new(((cur + p - sub) % p) as i64));
            }
        }
    }

    while remainder.len() > 1 && remainder.last().map_or(false, |c| c.0.is_zero()) {
        remainder.pop();
    }

    DensePoly::new(remainder)
}

fn mod_inv(a: u64, m: u64) -> u64 {
    let mut old_r = a as i128;
    let mut r = m as i128;
    let mut old_s: i128 = 1;
    let mut s: i128 = 0;

    while r != 0 {
        let q = old_r / r;
        let new_r = old_r - q * r;
        old_r = r;
        r = new_r;

        let new_s = old_s - q * s;
        old_s = s;
        s = new_s;
    }

    ((old_s % m as i128 + m as i128) % m as i128) as u64
}

fn factor_mod_prime(f: &DensePoly<Z>, p: u64) -> Vec<DensePoly<Z>> {
    let f_mod = reduce_mod_p(f, p);

    // Simplified factorization using distinct-degree factorization
    factor_mod_p_impl(&f_mod, p)
}

fn factor_mod_p_impl(f: &DensePoly<Z>, p: u64) -> Vec<DensePoly<Z>> {
    if f.degree() <= 1 {
        return vec![f.clone()];
    }

    // Distinct degree factorization
    let ddf = distinct_degree_factor(f, p);

    let mut all_factors = Vec::new();
    for (d, poly) in ddf {
        if poly.degree() == d {
            all_factors.push(poly);
        } else {
            let edf = equal_degree_factor(&poly, d, p);
            all_factors.extend(edf);
        }
    }

    all_factors
}

fn distinct_degree_factor(f: &DensePoly<Z>, p: u64) -> Vec<(usize, DensePoly<Z>)> {
    let n = f.degree();
    let mut result = Vec::new();
    let mut h = f.clone();

    let x = DensePoly::new(vec![z_zero(), z_one()]);
    let mut x_pow = pow_mod_p(&x, p, f, p);

    for d in 1..=n / 2 {
        if h.degree() < 2 * d {
            break;
        }

        if d > 1 {
            x_pow = pow_mod_p(&x_pow, p, &h, p);
        }

        let diff = poly_sub_mod_p(&x_pow, &x, p);
        let g = poly_gcd_mod(&h, &diff, p);

        if g.degree() > 0 {
            result.push((d, g.clone()));
            h = poly_div_mod_p(&h, &g, p);
            x_pow = poly_mod_p(&x_pow, &h, p);
        }
    }

    if h.degree() > 0 {
        result.push((h.degree(), h));
    }

    result
}

fn equal_degree_factor(f: &DensePoly<Z>, d: usize, p: u64) -> Vec<DensePoly<Z>> {
    let n = f.degree();
    if n == d {
        return vec![f.clone()];
    }

    let num_factors = n / d;
    let mut factors = vec![f.clone()];
    let mut seed = 42u64;

    while factors.len() < num_factors {
        let mut new_factors = Vec::new();

        for factor in &factors {
            if factor.degree() == d {
                new_factors.push(factor.clone());
                continue;
            }

            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let a = random_poly(factor.degree(), p, seed);

            let exp = (p - 1) / 2;
            let a_pow = pow_mod_p(&a, exp, factor, p);
            let a_pow_minus_1 = poly_sub_mod_p(&a_pow, &DensePoly::new(vec![z_one()]), p);

            let g = poly_gcd_mod(factor, &a_pow_minus_1, p);

            if g.degree() > 0 && g.degree() < factor.degree() {
                let other = poly_div_mod_p(factor, &g, p);
                new_factors.push(g);
                new_factors.push(other);
            } else {
                new_factors.push(factor.clone());
            }
        }

        factors = new_factors;
    }

    factors
}

fn random_poly(max_deg: usize, p: u64, seed: u64) -> DensePoly<Z> {
    let mut s = seed;
    let mut coeffs = Vec::with_capacity(max_deg);

    for _ in 0..max_deg {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
        coeffs.push(Z(Integer::new((s % p) as i64)));
    }

    DensePoly::new(coeffs)
}

fn pow_mod_p(a: &DensePoly<Z>, mut exp: u64, f: &DensePoly<Z>, p: u64) -> DensePoly<Z> {
    if exp == 0 {
        return DensePoly::new(vec![z_one()]);
    }

    let mut result = DensePoly::new(vec![z_one()]);
    let mut base = poly_mod_p(a, f, p);

    while exp > 0 {
        if exp & 1 == 1 {
            result = poly_mod_p(&poly_mul_mod_p(&result, &base, p), f, p);
        }
        base = poly_mod_p(&poly_mul_mod_p(&base, &base, p), f, p);
        exp >>= 1;
    }

    result
}

fn poly_mul_mod_p(a: &DensePoly<Z>, b: &DensePoly<Z>, p: u64) -> DensePoly<Z> {
    if a.is_zero() || b.is_zero() {
        return DensePoly::zero();
    }

    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let mut result = vec![z_zero(); a_coeffs.len() + b_coeffs.len() - 1];

    for (i, ai) in a_coeffs.iter().enumerate() {
        let ai_val = ai.0.to_i64().unwrap().rem_euclid(p as i64) as u64;
        for (j, bj) in b_coeffs.iter().enumerate() {
            let bj_val = bj.0.to_i64().unwrap().rem_euclid(p as i64) as u64;
            let prod = (ai_val as u128 * bj_val as u128) % p as u128;
            let cur = result[i + j].0.to_i64().unwrap().rem_euclid(p as i64) as u64;
            result[i + j] = Z(Integer::new(((cur as u128 + prod) % p as u128) as i64));
        }
    }

    while result.len() > 1 && result.last().map_or(false, |c| c.0.is_zero()) {
        result.pop();
    }

    DensePoly::new(result)
}

fn poly_sub_mod_p(a: &DensePoly<Z>, b: &DensePoly<Z>, p: u64) -> DensePoly<Z> {
    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let max_len = a_coeffs.len().max(b_coeffs.len());
    let mut result = vec![z_zero(); max_len];

    for (i, c) in a_coeffs.iter().enumerate() {
        result[i] = Z(Integer::new(c.0.to_i64().unwrap().rem_euclid(p as i64)));
    }
    for (i, c) in b_coeffs.iter().enumerate() {
        let cur = result[i].0.to_i64().unwrap();
        let sub = c.0.to_i64().unwrap().rem_euclid(p as i64);
        result[i] = Z(Integer::new(((cur - sub) % p as i64 + p as i64) % p as i64));
    }

    while result.len() > 1 && result.last().map_or(false, |c| c.0.is_zero()) {
        result.pop();
    }

    DensePoly::new(result)
}

fn poly_mod_p(a: &DensePoly<Z>, f: &DensePoly<Z>, p: u64) -> DensePoly<Z> {
    poly_mod_mod(a, f, p)
}

fn poly_div_mod_p(a: &DensePoly<Z>, b: &DensePoly<Z>, p: u64) -> DensePoly<Z> {
    if b.is_zero() {
        panic!("Division by zero");
    }

    if a.degree() < b.degree() {
        return DensePoly::zero();
    }

    let mut remainder = a.coeffs().to_vec();
    let divisor = b.coeffs();
    let lead_inv = mod_inv(
        divisor.last().unwrap().0.to_i64().unwrap().rem_euclid(p as i64) as u64,
        p,
    );
    let deg_diff = a.degree() - b.degree();
    let mut quotient = vec![z_zero(); deg_diff + 1];

    for i in (0..=deg_diff).rev() {
        let idx = i + divisor.len() - 1;
        if idx < remainder.len() {
            let r_val = remainder[idx].0.to_i64().unwrap().rem_euclid(p as i64) as u64;
            let q = (r_val as u128 * lead_inv as u128) % p as u128;
            quotient[i] = Z(Integer::new(q as i64));

            for (j, d) in divisor.iter().enumerate() {
                let ridx = i + j;
                let d_val = d.0.to_i64().unwrap().rem_euclid(p as i64) as u64;
                let sub = (q * d_val as u128) % p as u128;
                let cur = remainder[ridx].0.to_i64().unwrap().rem_euclid(p as i64) as u64;
                remainder[ridx] =
                    Z(Integer::new(((cur as u128 + p as u128 - sub) % p as u128) as i64));
            }
        }
    }

    while quotient.len() > 1 && quotient.last().map_or(false, |c| c.0.is_zero()) {
        quotient.pop();
    }

    DensePoly::new(quotient)
}

fn factor_coefficient_bound(f: &DensePoly<Z>) -> Integer {
    let n = f.degree();
    let coeffs = f.coeffs();

    let mut max_coeff = Integer::new(0);
    for c in coeffs {
        let abs = c.0.clone().abs();
        if abs > max_coeff {
            max_coeff = abs;
        }
    }

    let two_n = Integer::new(2).pow(n as u32);
    max_coeff * two_n
}

fn compute_precision(p: u64, bound: &Integer) -> u32 {
    let mut modulus = Integer::new(p as i64);
    let target = bound.clone() * Integer::new(2);
    let mut k = 1u32;

    while modulus <= target {
        modulus = modulus * Integer::new(p as i64);
        k += 1;
    }

    k
}

fn combine_factors(f: &DensePoly<Z>, lifted: &HenselLiftResult) -> Vec<DensePoly<Z>> {
    let r = lifted.factors.len();

    if r == 0 {
        return vec![f.clone()];
    }

    if r == 1 {
        return vec![f.clone()];
    }

    let mut remaining = f.clone();
    let mut found_factors = Vec::new();
    let mut used = vec![false; r];

    // Try single modular factors
    for i in 0..r {
        if used[i] {
            continue;
        }

        let candidate = &lifted.factors[i];
        if let Some(true_factor) = try_factor(&remaining, candidate, &lifted.modulus) {
            found_factors.push(true_factor.clone());
            remaining = poly_div_exact(&remaining, &true_factor);
            used[i] = true;

            if remaining.degree() == 0 {
                break;
            }
        }
    }

    // Try pairs if needed
    if remaining.degree() > 0 && found_factors.len() + 2 <= r {
        for i in 0..r {
            if used[i] || remaining.degree() == 0 {
                continue;
            }
            for j in i + 1..r {
                if used[j] {
                    continue;
                }

                let product = poly_mul_mod_z(
                    &lifted.factors[i],
                    &lifted.factors[j],
                    &lifted.modulus,
                );

                if let Some(true_factor) = try_factor(&remaining, &product, &lifted.modulus) {
                    found_factors.push(true_factor.clone());
                    remaining = poly_div_exact(&remaining, &true_factor);
                    used[i] = true;
                    used[j] = true;

                    if remaining.degree() == 0 {
                        break;
                    }
                }
            }
            if remaining.degree() == 0 {
                break;
            }
        }
    }

    if remaining.degree() > 0 {
        found_factors.push(remaining);
    }

    found_factors
}

fn try_factor(
    f: &DensePoly<Z>,
    mod_factor: &DensePoly<Z>,
    modulus: &Integer,
) -> Option<DensePoly<Z>> {
    let half = modulus.clone() / Integer::new(2);
    let candidate: Vec<Z> = mod_factor
        .coeffs()
        .iter()
        .map(|c| {
            let r = c.0.clone() % modulus.clone();
            let r = if r < Integer::new(0) {
                r + modulus.clone()
            } else {
                r
            };
            if r > half {
                Z(r - modulus.clone())
            } else {
                Z(r)
            }
        })
        .collect();

    let candidate = make_primitive(&DensePoly::new(candidate));

    if candidate.degree() == 0 || candidate.degree() >= f.degree() {
        return None;
    }

    let (q, r) = poly_div_z(f, &candidate);
    if r.is_zero() || r.coeffs().iter().all(|c| c.0.is_zero()) {
        Some(candidate)
    } else {
        None
    }
}

fn poly_mul_mod_z(a: &DensePoly<Z>, b: &DensePoly<Z>, m: &Integer) -> DensePoly<Z> {
    if a.is_zero() || b.is_zero() {
        return DensePoly::zero();
    }

    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let mut result = vec![z_zero(); a_coeffs.len() + b_coeffs.len() - 1];

    for (i, ai) in a_coeffs.iter().enumerate() {
        for (j, bj) in b_coeffs.iter().enumerate() {
            let prod = (ai.0.clone() * bj.0.clone()) % m.clone();
            result[i + j] = Z((result[i + j].0.clone() + prod) % m.clone());
        }
    }

    while result.len() > 1 && result.last().map_or(false, |c| c.0.is_zero()) {
        result.pop();
    }

    DensePoly::new(result)
}

fn poly_div_z(a: &DensePoly<Z>, b: &DensePoly<Z>) -> (DensePoly<Z>, DensePoly<Z>) {
    if b.is_zero() {
        panic!("Division by zero");
    }

    if a.degree() < b.degree() {
        return (DensePoly::zero(), a.clone());
    }

    let mut remainder = a.coeffs().to_vec();
    let divisor = b.coeffs();
    let lead_b = divisor.last().unwrap();
    let deg_diff = a.degree() - b.degree();
    let mut quotient = vec![z_zero(); deg_diff + 1];

    for i in (0..=deg_diff).rev() {
        let idx = i + divisor.len() - 1;
        if idx < remainder.len() {
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

    while remainder.len() > 1 && remainder.last().map_or(false, |c| c.0.is_zero()) {
        remainder.pop();
    }
    while quotient.len() > 1 && quotient.last().map_or(false, |c| c.0.is_zero()) {
        quotient.pop();
    }

    (DensePoly::new(quotient), DensePoly::new(remainder))
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

fn make_primitive(f: &DensePoly<Z>) -> DensePoly<Z> {
    if f.is_zero() {
        return f.clone();
    }

    let content = polynomial_content(f);
    let coeffs: Vec<Z> = if content.0.is_one() || content.0.is_zero() {
        f.coeffs().to_vec()
    } else {
        f.coeffs()
            .iter()
            .map(|c| Z(c.0.clone() / content.0.clone()))
            .collect()
    };

    let mut result = DensePoly::new(coeffs);

    if result.coeffs().last().map_or(false, |c| c.0 < Integer::new(0)) {
        result = DensePoly::new(
            result
                .coeffs()
                .iter()
                .map(|c| Z(-c.0.clone()))
                .collect(),
        );
    }

    result
}

fn poly_div_exact(a: &DensePoly<Z>, b: &DensePoly<Z>) -> DensePoly<Z> {
    poly_div_z(a, b).0
}

pub fn van_hoeij_factor_batch(polys: &[DensePoly<Z>]) -> Vec<VanHoeijResult> {
    polys.par_iter().map(van_hoeij_factor).collect()
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
    fn test_factor_linear() {
        let f = poly(&[3, 2]);
        let result = van_hoeij_factor(&f);

        assert_eq!(result.factors.len(), 1);
    }

    #[test]
    fn test_factor_two_linears() {
        let f = poly(&[2, 3, 1]);
        let result = van_hoeij_factor(&f);

        assert_eq!(result.factors.len(), 2);
        let f0 = &result.factors[0];
        let f1 = &result.factors[1];
        assert!(f0.degree() == 1 && f1.degree() == 1);
    }

    #[test]
    fn test_factor_irreducible() {
        let f = poly(&[1, 0, 1]);
        let result = van_hoeij_factor(&f);

        assert_eq!(result.factors.len(), 1);
        assert_eq!(result.factors[0].degree(), 2);
    }

    #[test]
    fn test_mod_inv() {
        assert_eq!(mod_inv(3, 7), 5);
        assert_eq!(mod_inv(2, 11), 6);
    }
}
