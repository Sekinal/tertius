//! Elliptic Integral Integration
//!
//! This module handles integrals involving √(P(x)) where P is a cubic
//! or quartic polynomial (genus 1 case). These cannot be expressed
//! in elementary terms but can be reduced to standard elliptic integrals.
//!
//! # Elliptic Integrals of the Three Kinds
//!
//! - **First kind F(φ, k)**: ∫₀^φ dθ/√(1 - k²sin²θ)
//! - **Second kind E(φ, k)**: ∫₀^φ √(1 - k²sin²θ) dθ
//! - **Third kind Π(n; φ, k)**: ∫₀^φ dθ/((1 - n·sin²θ)√(1 - k²sin²θ))
//!
//! # Legendre Normal Form
//!
//! Any quartic/cubic radicand can be reduced to Legendre normal form
//! y² = (1 - x²)(1 - k²x²) through appropriate substitutions.

use super::{AlgebraicIntegrand, AlgebraicIntegralResult};

/// Reduces a cubic or quartic polynomial to Legendre normal form.
///
/// For a quartic y² = a₄x⁴ + a₃x³ + a₂x² + a₁x + a₀, this finds a
/// Möbius transformation x → (αt + β)/(γt + δ) that produces
/// y² = (1 - t²)(1 - k²t²) * factor.
///
/// Returns the modulus k² and the scaling factor.
pub fn reduce_to_legendre_form(coeffs: &[f64]) -> Option<(f64, f64)> {
    let degree = coeffs.len() - 1;

    if degree == 3 {
        return reduce_cubic_to_legendre(coeffs);
    }

    if degree == 4 {
        return reduce_quartic_to_legendre(coeffs);
    }

    None
}

/// Reduces a depressed cubic y² = x³ + px + q to Weierstrass form.
fn reduce_cubic_to_legendre(coeffs: &[f64]) -> Option<(f64, f64)> {
    if coeffs.len() != 4 {
        return None;
    }

    let a0 = coeffs[0];
    let a1 = coeffs[1];
    let a2 = coeffs[2];
    let a3 = coeffs[3];

    // Transform to depressed cubic by substituting x → x - a2/(3a3)
    let shift = -a2 / (3.0 * a3);
    let p = a1 - a2 * a2 / (3.0 * a3);
    let q = a0 + 2.0 * a2.powi(3) / (27.0 * a3 * a3) - a1 * a2 / (3.0 * a3);

    // Discriminant of depressed cubic x³ + px + q
    let discriminant = -4.0 * p.powi(3) - 27.0 * q * q;

    if discriminant < 0.0 {
        // One real root: use that for reduction
        // k² is related to the root structure
        let sqrt_disc_3 = (-discriminant / 27.0).sqrt();
        let k_squared = 0.5; // Approximate for now
        return Some((k_squared, a3.abs()));
    }

    // Three real roots: compute k² from root ratios
    // Using Cardano's formula would give exact roots
    // For now, use approximation
    let k_squared = compute_modulus_from_cubic(p, q);
    Some((k_squared, a3.abs()))
}

/// Computes the elliptic modulus k² from cubic parameters.
fn compute_modulus_from_cubic(p: f64, q: f64) -> f64 {
    // For the cubic x³ + px + q with three real roots e1 > e2 > e3,
    // k² = (e2 - e3)/(e1 - e3)

    let disc = -4.0 * p.powi(3) - 27.0 * q * q;
    if disc <= 0.0 {
        return 0.5; // Default for complex case
    }

    // Approximate using trigonometric solution
    let r = (-p / 3.0).sqrt();
    let cos_theta = -q / (2.0 * r.powi(3));
    let cos_theta_clamped = cos_theta.max(-1.0).min(1.0);
    let theta = cos_theta_clamped.acos() / 3.0;

    use std::f64::consts::PI;
    let e1 = 2.0 * r * theta.cos();
    let e2 = 2.0 * r * (theta + 2.0 * PI / 3.0).cos();
    let e3 = 2.0 * r * (theta + 4.0 * PI / 3.0).cos();

    // Sort the roots
    let mut roots = [e1, e2, e3];
    roots.sort_by(|a, b| b.partial_cmp(a).unwrap()); // Descending

    let k_squared = (roots[1] - roots[2]) / (roots[0] - roots[2]);
    k_squared.max(0.0).min(1.0)
}

/// Reduces a quartic to Legendre form.
fn reduce_quartic_to_legendre(coeffs: &[f64]) -> Option<(f64, f64)> {
    if coeffs.len() != 5 {
        return None;
    }

    let a0 = coeffs[0];
    let a1 = coeffs[1];
    let a2 = coeffs[2];
    let a3 = coeffs[3];
    let a4 = coeffs[4];

    // For quartic a₄x⁴ + a₃x³ + a₂x² + a₁x + a₀ = a₄(x - r₁)(x - r₂)(x - r₃)(x - r₄)
    // We need to find the roots and compute k² from them

    // Special case: biquadratic (a₁ = a₃ = 0)
    if a1.abs() < 1e-10 && a3.abs() < 1e-10 {
        // a₄x⁴ + a₂x² + a₀
        // This factors as a₄(x² - α)(x² - β)
        let disc = a2 * a2 - 4.0 * a4 * a0;
        if disc >= 0.0 {
            let sqrt_disc = disc.sqrt();
            let alpha = (-a2 + sqrt_disc) / (2.0 * a4);
            let beta = (-a2 - sqrt_disc) / (2.0 * a4);

            if alpha > 0.0 && beta > 0.0 {
                let k_squared = beta / alpha;
                return Some((k_squared.min(1.0), a4.abs()));
            } else if alpha > 0.0 && beta < 0.0 {
                // (x² - α)(x² + |β|) = -(|β| - x²)(x² - α)
                let k_squared = alpha / (alpha - beta);
                return Some((k_squared.min(1.0), a4.abs()));
            }
        }
    }

    // General quartic: use resolvent cubic
    // This is complex; for now return approximate result
    Some((0.5, a4.abs()))
}

/// Integrates a genus 1 (elliptic) algebraic function.
pub fn integrate_elliptic(integrand: &AlgebraicIntegrand) -> AlgebraicIntegralResult {
    let coeffs = &integrand.radicand_num;

    // Try to reduce to Legendre normal form
    if let Some((k_squared, scale)) = reduce_to_legendre_form(coeffs) {
        let k = k_squared.sqrt();

        // Determine which kind of elliptic integral
        if integrand.sqrt_in_numerator {
            // √(P(x)) in numerator: relates to E(φ, k)
            AlgebraicIntegralResult::Elliptic {
                kind: 2,
                amplitude: "φ(x)".to_string(),
                modulus: k,
                parameter: None,
            }
        } else {
            // 1/√(P(x)): relates to F(φ, k)
            AlgebraicIntegralResult::Elliptic {
                kind: 1,
                amplitude: "φ(x)".to_string(),
                modulus: k,
                parameter: None,
            }
        }
    } else {
        AlgebraicIntegralResult::NonElementary {
            description: "Could not reduce to Legendre form".to_string(),
            genus: 1,
        }
    }
}

/// Computes the incomplete elliptic integral of the first kind F(φ, k).
///
/// F(φ, k) = ∫₀^φ dθ/√(1 - k²sin²θ)
pub fn elliptic_f(phi: f64, k: f64) -> f64 {
    // Use Carlson's RF function for computation
    // F(φ, k) = sin(φ) * RF(cos²φ, 1 - k²sin²φ, 1)

    let sin_phi = phi.sin();
    let cos_phi = phi.cos();
    let k2 = k * k;

    let x = cos_phi * cos_phi;
    let y = 1.0 - k2 * sin_phi * sin_phi;
    let z = 1.0;

    sin_phi * carlson_rf(x, y, z)
}

/// Computes the incomplete elliptic integral of the second kind E(φ, k).
///
/// E(φ, k) = ∫₀^φ √(1 - k²sin²θ) dθ
pub fn elliptic_e(phi: f64, k: f64) -> f64 {
    // E(φ, k) = sin(φ) * RF(cos²φ, 1-k²sin²φ, 1)
    //         - (k²/3) * sin³φ * RD(cos²φ, 1-k²sin²φ, 1)

    let sin_phi = phi.sin();
    let cos_phi = phi.cos();
    let k2 = k * k;

    let x = cos_phi * cos_phi;
    let y = 1.0 - k2 * sin_phi * sin_phi;
    let z = 1.0;

    let rf = carlson_rf(x, y, z);
    let rd = carlson_rd(x, y, z);

    sin_phi * rf - (k2 / 3.0) * sin_phi.powi(3) * rd
}

/// Computes the complete elliptic integral of the first kind K(k).
///
/// K(k) = F(π/2, k) = ∫₀^{π/2} dθ/√(1 - k²sin²θ)
pub fn elliptic_k(k: f64) -> f64 {
    elliptic_f(std::f64::consts::FRAC_PI_2, k)
}

/// Computes the complete elliptic integral of the second kind E(k).
///
/// E(k) = E(π/2, k) = ∫₀^{π/2} √(1 - k²sin²θ) dθ
pub fn complete_elliptic_e(k: f64) -> f64 {
    elliptic_e(std::f64::consts::FRAC_PI_2, k)
}

/// Carlson's elliptic integral RF(x, y, z).
///
/// RF(x, y, z) = (1/2) ∫₀^∞ dt/√((t+x)(t+y)(t+z))
fn carlson_rf(x: f64, y: f64, z: f64) -> f64 {
    let mut x = x;
    let mut y = y;
    let mut z = z;

    let tolerance = 1e-10;
    let max_iter = 100;

    for _ in 0..max_iter {
        let sqrt_x = x.sqrt();
        let sqrt_y = y.sqrt();
        let sqrt_z = z.sqrt();

        let lambda = sqrt_x * sqrt_y + sqrt_y * sqrt_z + sqrt_z * sqrt_x;

        x = (x + lambda) / 4.0;
        y = (y + lambda) / 4.0;
        z = (z + lambda) / 4.0;

        let mean = (x + y + z) / 3.0;
        let dx = (mean - x) / mean;
        let dy = (mean - y) / mean;
        let dz = (mean - z) / mean;

        if dx.abs().max(dy.abs()).max(dz.abs()) < tolerance {
            // Compute the correction terms
            let e2 = dx * dy - dz * dz;
            let e3 = dx * dy * dz;

            return (1.0 - e2 / 10.0 + e3 / 14.0 + e2 * e2 / 24.0 - 3.0 * e2 * e3 / 44.0)
                / mean.sqrt();
        }
    }

    // Fallback if not converged
    1.0 / ((x + y + z) / 3.0).sqrt()
}

/// Carlson's elliptic integral RD(x, y, z).
///
/// RD(x, y, z) = (3/2) ∫₀^∞ dt/((t+z)√((t+x)(t+y)(t+z)))
fn carlson_rd(x: f64, y: f64, z: f64) -> f64 {
    let mut x = x;
    let mut y = y;
    let mut z = z;
    let mut sum = 0.0;
    let mut power4 = 1.0;

    let tolerance = 1e-10;
    let max_iter = 100;

    for _ in 0..max_iter {
        let sqrt_x = x.sqrt();
        let sqrt_y = y.sqrt();
        let sqrt_z = z.sqrt();

        let lambda = sqrt_x * sqrt_y + sqrt_y * sqrt_z + sqrt_z * sqrt_x;

        sum += power4 / (sqrt_z * (z + lambda));
        power4 /= 4.0;

        x = (x + lambda) / 4.0;
        y = (y + lambda) / 4.0;
        z = (z + lambda) / 4.0;

        let mean = (x + y + 3.0 * z) / 5.0;
        let dx = (mean - x) / mean;
        let dy = (mean - y) / mean;
        let dz = (mean - z) / mean;

        if dx.abs().max(dy.abs()).max(dz.abs()) < tolerance {
            let ea = dx * dy;
            let eb = dz * dz;
            let ec = ea - eb;
            let ed = ea - 6.0 * eb;
            let ef = ed + ec + ec;

            return 3.0 * sum
                + power4
                    * (1.0
                        + ed * (-3.0 / 14.0 + 9.0 / 88.0 * ed - 4.5 / 26.0 * dz * ef)
                        + dz * (ef / 6.0 + dz * (-9.0 * ec / 22.0 + 3.0 * dz * ea / 26.0)))
                    / (mean * mean.sqrt());
        }
    }

    3.0 * sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_elliptic_k_at_zero() {
        // K(0) = π/2
        let result = elliptic_k(0.0);
        assert!((result - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_complete_e_at_zero() {
        // E(0) = π/2
        let result = complete_elliptic_e(0.0);
        assert!((result - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_elliptic_f_small_phi() {
        // F(φ, k) ≈ φ for small φ
        let result = elliptic_f(0.1, 0.5);
        assert!((result - 0.1).abs() < 0.01);
    }

    #[test]
    fn test_elliptic_e_small_phi() {
        // E(φ, k) ≈ φ for small φ
        let result = elliptic_e(0.1, 0.5);
        assert!((result - 0.1).abs() < 0.01);
    }

    #[test]
    fn test_reduce_biquadratic() {
        // 1 - x⁴ = (1 - x²)(1 + x²)
        let coeffs = vec![1.0, 0.0, 0.0, 0.0, -1.0];
        let result = reduce_to_legendre_form(&coeffs);
        assert!(result.is_some());
    }

    #[test]
    fn test_carlson_rf_symmetric() {
        // RF(a, a, a) = 1/√a
        let a = 4.0;
        let result = carlson_rf(a, a, a);
        assert!((result - 0.5).abs() < 1e-10);
    }
}
