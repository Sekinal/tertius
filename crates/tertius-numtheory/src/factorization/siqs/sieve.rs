//! Sieving for SIQS.
//!
//! The sieve finds values x in [-M, M] where Q(x) is smooth over the factor base.

use tertius_integers::Integer;

use super::{Relation, SiqsParams};

/// Sieve for smooth relations.
pub fn sieve_for_relations(n: &Integer, factor_base: &[u64], params: &SiqsParams) -> Vec<Relation> {
    let mut relations = Vec::new();

    // Simple implementation: just try values and check for smoothness
    // A real implementation would use proper sieving

    let sqrt_n = integer_sqrt(n);
    let m = params.sieve_interval as i64;

    for offset in 0..100 {
        // Try x = sqrt(n) + offset
        let x = sqrt_n.clone() + Integer::from(offset);
        let x_sq = x.clone() * x.clone();
        let target = x_sq - n.clone();

        if let Some((exponents, factorization)) = try_factor_smooth(&target, factor_base) {
            relations.push(Relation {
                x,
                exponents,
                factorization,
            });

            // Check if we have enough relations
            if relations.len() > factor_base.len() + 10 {
                return relations;
            }
        }

        // Also try negative offset
        if offset > 0 {
            let x = sqrt_n.clone() - Integer::from(offset);
            if x > Integer::from(1i64) {
                let x_sq = x.clone() * x.clone();
                let target = if x_sq >= n.clone() {
                    x_sq - n.clone()
                } else {
                    n.clone() - x_sq
                };

                if let Some((exponents, factorization)) = try_factor_smooth(&target, factor_base) {
                    relations.push(Relation {
                        x,
                        exponents,
                        factorization,
                    });

                    if relations.len() > factor_base.len() + 10 {
                        return relations;
                    }
                }
            }
        }
    }

    relations
}

/// Try to factor target over the factor base.
/// Returns the exponent vector (mod 2) and full factorization if smooth.
fn try_factor_smooth(
    target: &Integer,
    factor_base: &[u64],
) -> Option<(Vec<u8>, Vec<(usize, u32)>)> {
    use num_traits::{One, Zero};

    if target.is_zero() {
        return None;
    }

    let mut remaining = target.abs();
    let mut exponents = Vec::new();
    let mut factorization = Vec::new();

    for (idx, &p) in factor_base.iter().enumerate() {
        let p_int = Integer::from(p);
        let mut exp = 0u32;

        while (remaining.clone() % p_int.clone()).is_zero() {
            remaining = remaining / p_int.clone();
            exp += 1;
        }

        if exp > 0 {
            exponents.push((exp % 2) as u8);
            factorization.push((idx, exp));
        }

        if remaining.is_one() {
            return Some((exponents, factorization));
        }
    }

    // Check if remaining is 1 (fully factored)
    if remaining.is_one() {
        Some((exponents, factorization))
    } else {
        None
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
