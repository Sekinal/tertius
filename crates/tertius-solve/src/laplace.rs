//! Laplace transform starter utilities.

/// Supported time-domain basis terms.
#[derive(Clone, Debug, PartialEq)]
pub enum TimeTerm {
    /// c
    Constant(f64),
    /// t^n
    Power(u32),
    /// e^{a t}
    Exp { a: f64 },
    /// sin(omega t)
    Sin { omega: f64 },
    /// cos(omega t)
    Cos { omega: f64 },
}

/// Rational image in Laplace domain: numerator/denominator in ascending powers of s.
#[derive(Clone, Debug, PartialEq)]
pub struct LaplaceImage {
    /// Numerator coefficients.
    pub numerator: Vec<f64>,
    /// Denominator coefficients.
    pub denominator: Vec<f64>,
}

impl LaplaceImage {
    /// Evaluate F(s).
    #[must_use]
    pub fn eval(&self, s: f64) -> f64 {
        let num = eval_poly(&self.numerator, s);
        let den = eval_poly(&self.denominator, s);
        num / den
    }
}

/// Laplace transform for supported basis terms.
#[must_use]
pub fn laplace_transform_basic(term: &TimeTerm) -> LaplaceImage {
    match term {
        TimeTerm::Constant(c) => LaplaceImage {
            numerator: vec![*c],
            denominator: vec![0.0, 1.0],
        },
        TimeTerm::Power(n) => {
            let mut den = vec![0.0; (*n as usize) + 2];
            den[(*n as usize) + 1] = 1.0;
            LaplaceImage {
                numerator: vec![factorial(*n)],
                denominator: den,
            }
        }
        TimeTerm::Exp { a } => LaplaceImage {
            numerator: vec![1.0],
            denominator: vec![-*a, 1.0],
        },
        TimeTerm::Sin { omega } => LaplaceImage {
            numerator: vec![*omega],
            denominator: vec![omega * omega, 0.0, 1.0],
        },
        TimeTerm::Cos { omega } => LaplaceImage {
            numerator: vec![0.0, 1.0],
            denominator: vec![omega * omega, 0.0, 1.0],
        },
    }
}

/// Inverse Laplace transform for supported canonical images.
#[must_use]
pub fn inverse_laplace_basic(image: &LaplaceImage) -> Option<TimeTerm> {
    if image.numerator.len() == 1 && image.denominator == [0.0, 1.0] {
        return Some(TimeTerm::Constant(image.numerator[0]));
    }
    if image.numerator.len() == 1
        && image.denominator.len() == 2
        && approx_eq(image.denominator[1], 1.0)
        && !approx_eq(image.denominator[0], 0.0)
        && approx_eq(image.numerator[0], 1.0)
    {
        return Some(TimeTerm::Exp {
            a: -image.denominator[0],
        });
    }
    if image.numerator.len() == 1
        && image.denominator.len() == 3
        && approx_eq(image.denominator[1], 0.0)
        && approx_eq(image.denominator[2], 1.0)
        && approx_eq(
            image.denominator[0],
            image.numerator[0] * image.numerator[0],
        )
    {
        return Some(TimeTerm::Sin {
            omega: image.numerator[0],
        });
    }
    if image.numerator == [0.0, 1.0]
        && image.denominator.len() == 3
        && approx_eq(image.denominator[1], 0.0)
        && approx_eq(image.denominator[2], 1.0)
    {
        return Some(TimeTerm::Cos {
            omega: image.denominator[0].sqrt(),
        });
    }
    None
}

fn eval_poly(coeffs: &[f64], s: f64) -> f64 {
    coeffs
        .iter()
        .enumerate()
        .map(|(i, c)| c * s.powi(i as i32))
        .sum()
}

fn factorial(n: u32) -> f64 {
    (1..=n).fold(1.0, |acc, k| acc * f64::from(k))
}

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

/// Laplace-domain image `Y(s)` for `y' + a y = b`, `y(0)=y0`.
///
/// Formula:
/// `Y(s) = (y0*s + a*y0 + b) / (s^2 + a*s)`.
#[must_use]
pub fn first_order_linear_ivp_image(a: f64, b: f64, y0: f64) -> LaplaceImage {
    LaplaceImage {
        // y0*s + a*y0 + b
        numerator: vec![a * y0 + b, y0],
        // s^2 + a*s
        denominator: vec![0.0, a, 1.0],
    }
}

/// Laplace-domain image `Y(s)` for `y'' + b y' + c y = 0`,
/// with `y(0)=y0`, `y'(0)=v0`.
///
/// Formula:
/// `Y(s) = (s*y0 + v0 + b*y0) / (s^2 + b*s + c)`.
#[must_use]
pub fn second_order_homogeneous_ivp_image(b: f64, c: f64, y0: f64, v0: f64) -> LaplaceImage {
    LaplaceImage {
        // s*y0 + v0 + b*y0
        numerator: vec![v0 + b * y0, y0],
        // s^2 + b*s + c
        denominator: vec![c, b, 1.0],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_laplace_sin() {
        let image = laplace_transform_basic(&TimeTerm::Sin { omega: 2.0 });
        // L{sin(2t)} = 2/(s^2 + 4)
        let expected = 2.0 / (3.0f64.powi(2) + 4.0);
        assert!((image.eval(3.0) - expected).abs() < 1e-9);
    }

    #[test]
    fn test_laplace_exp() {
        let image = laplace_transform_basic(&TimeTerm::Exp { a: 2.0 });
        // L{e^{2t}} = 1/(s-2)
        assert!((image.eval(5.0) - (1.0 / 3.0)).abs() < 1e-9);
    }

    #[test]
    fn test_inverse_laplace_sin() {
        let img = LaplaceImage {
            numerator: vec![2.0],
            denominator: vec![4.0, 0.0, 1.0],
        };
        assert_eq!(
            inverse_laplace_basic(&img),
            Some(TimeTerm::Sin { omega: 2.0 })
        );
    }

    #[test]
    fn test_first_order_linear_ivp_image_shape() {
        // y' + 2y = 4, y(0)=1 => Y(s) = (s+6)/(s^2+2s)
        let image = first_order_linear_ivp_image(2.0, 4.0, 1.0);
        assert_eq!(image.numerator, vec![6.0, 1.0]);
        assert_eq!(image.denominator, vec![0.0, 2.0, 1.0]);
    }

    #[test]
    fn test_second_order_homogeneous_ivp_image_shape() {
        // y'' + y = 0, y(0)=0, y'(0)=1 => Y(s)=1/(s^2+1)
        let image = second_order_homogeneous_ivp_image(0.0, 1.0, 0.0, 1.0);
        assert_eq!(image.numerator, vec![1.0, 0.0]);
        assert_eq!(image.denominator, vec![1.0, 0.0, 1.0]);
    }
}
