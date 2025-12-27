//! F5 criteria for detecting useless S-polynomials.
//!
//! These criteria allow early rejection of pairs that would reduce to zero,
//! significantly improving performance.

use crate::labeled_poly::LabeledPoly;
use crate::monomial::PackedMonomial;
use crate::signature::{SignedPair, Signature};

/// Checks the F5 criterion: whether a signature is rewritable.
///
/// A signature σ is rewritable if there exists a polynomial g in the basis
/// such that lm(g) divides σ.monomial and sig(g) < σ with the same index.
pub fn is_rewritable<R: Clone + PartialEq>(
    sig: &Signature,
    basis: &[LabeledPoly<R>],
) -> bool {
    for g in basis {
        if g.signature.index != sig.index {
            continue;
        }

        if let Some(lm) = g.leading_monomial() {
            // Check if lm(g) divides the monomial part of σ
            if sig.monomial.is_divisible_by(lm) {
                // Check if sig(g) < σ
                if g.signature.cmp_pot(sig) == std::cmp::Ordering::Less {
                    return true;
                }
            }
        }
    }
    false
}

/// Checks Buchberger's first criterion (product criterion).
///
/// If lcm(lm(f), lm(g)) = lm(f) * lm(g) (i.e., leading monomials are coprime),
/// then S(f, g) reduces to zero.
pub fn product_criterion(lm_f: &PackedMonomial, lm_g: &PackedMonomial) -> bool {
    lm_f.is_coprime(lm_g)
}

/// Checks the chain criterion (Buchberger's second criterion).
///
/// If there exists h in the basis such that:
/// - lm(h) divides lcm(lm(f), lm(g))
/// - The pairs (f, h) and (g, h) have already been processed
/// then the pair (f, g) is redundant.
pub fn chain_criterion<R: Clone + PartialEq>(
    pair: &SignedPair,
    basis: &[LabeledPoly<R>],
    processed_pairs: &[(usize, usize)],
) -> bool {
    for (k, h) in basis.iter().enumerate() {
        if k == pair.i || k == pair.j {
            continue;
        }

        if let Some(lm_h) = h.leading_monomial() {
            // Check if lm(h) divides lcm
            if pair.lcm.is_divisible_by(lm_h) {
                // Check if pairs (i, k) and (j, k) are already processed
                let pair_ik = if pair.i < k { (pair.i, k) } else { (k, pair.i) };
                let pair_jk = if pair.j < k { (pair.j, k) } else { (k, pair.j) };

                if processed_pairs.contains(&pair_ik) && processed_pairs.contains(&pair_jk) {
                    return true;
                }
            }
        }
    }
    false
}

/// Checks the syzygy criterion (F5 criterion).
///
/// If the signature corresponds to a known syzygy (relation among generators),
/// then the S-polynomial reduces to zero.
///
/// In F5, principal syzygies are of the form: g_i * lm(g_j) - g_j * lm(g_i) = 0
pub fn syzygy_criterion<R: Clone + PartialEq>(
    pair: &SignedPair,
    basis: &[LabeledPoly<R>],
) -> bool {
    // For pairs where i < j, check if the signature indicates a principal syzygy
    if pair.i >= pair.j {
        return false;
    }

    // Get leading monomials
    let lm_i = basis.get(pair.i).and_then(|p| p.leading_monomial());
    let lm_j = basis.get(pair.j).and_then(|p| p.leading_monomial());

    if let (Some(lm_i), Some(lm_j)) = (lm_i, lm_j) {
        // Check if this is a principal syzygy
        // The signature of a principal syzygy g_i * lm(g_j) has index i
        // and monomial lm(g_j)
        if pair.signature.index == pair.i {
            // The signature monomial should be divisible by lm(g_j)
            if pair.signature.monomial.is_divisible_by(lm_j) {
                // Check the ratio
                if let Some(ratio) = pair.signature.monomial.div(lm_j) {
                    // If ratio is 1, this is exactly the principal syzygy
                    if ratio.is_one() {
                        return true;
                    }
                }
            }
        }
    }

    false
}

/// Selects pairs to process using the "sugar" strategy.
///
/// Returns pairs sorted by increasing sugar degree, which tends to
/// produce smaller intermediate polynomials.
pub fn sugar_selection(pairs: &mut [SignedPair]) -> Vec<SignedPair> {
    // Sort by sugar, then by signature
    pairs.sort_by(|a, b| {
        a.sugar.cmp(&b.sugar).then_with(|| a.signature.cmp_pot(&b.signature))
    });

    // Take all pairs with minimum sugar
    if pairs.is_empty() {
        return vec![];
    }

    let min_sugar = pairs[0].sugar;
    pairs
        .iter()
        .take_while(|p| p.sugar == min_sugar)
        .cloned()
        .collect()
}

/// Selects pairs by signature (F5 strategy).
///
/// Returns pairs with the smallest signatures first.
pub fn signature_selection(pairs: &mut [SignedPair]) -> Vec<SignedPair> {
    pairs.sort_by(|a, b| a.signature.cmp_pot(&b.signature));

    if pairs.is_empty() {
        return vec![];
    }

    // Take pairs with the smallest signature
    let min_sig = pairs[0].signature.clone();
    pairs
        .iter()
        .take_while(|p| p.signature.cmp_pot(&min_sig) == std::cmp::Ordering::Equal)
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_product_criterion() {
        let m1 = PackedMonomial::new(&[2, 0, 0]); // x^2
        let m2 = PackedMonomial::new(&[0, 3, 0]); // y^3

        assert!(product_criterion(&m1, &m2)); // Coprime

        let m3 = PackedMonomial::new(&[1, 1, 0]); // xy
        assert!(!product_criterion(&m1, &m3)); // Not coprime (share x)
    }
}
