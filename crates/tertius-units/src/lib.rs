//! Units and dimensional analysis primitives.

/// A physical dimension represented in SI base exponents.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Dimension {
    /// Exponents for [L, M, T, I, Theta, N, J].
    pub exponents: [i8; 7],
}

impl Dimension {
    /// Returns the dimensionless dimension.
    #[must_use]
    pub fn dimensionless() -> Self {
        Self { exponents: [0; 7] }
    }

    /// Raises this dimension to an integer power.
    #[must_use]
    pub fn powi(self, power: i8) -> Self {
        let mut exponents = self.exponents;
        for exp in &mut exponents {
            *exp *= power;
        }
        Self { exponents }
    }
}

/// A concrete unit with a scale relative to SI base units.
#[derive(Clone, Debug, PartialEq)]
pub struct Unit {
    /// Human-readable symbol.
    pub symbol: &'static str,
    /// Scale factor relative to SI base units.
    pub scale_to_si: f64,
    /// Underlying physical dimension.
    pub dimension: Dimension,
}

impl Unit {
    /// Creates a custom unit.
    #[must_use]
    pub fn new(symbol: &'static str, scale_to_si: f64, dimension: Dimension) -> Self {
        Self {
            symbol,
            scale_to_si,
            dimension,
        }
    }

    /// Meter (m).
    #[must_use]
    pub fn meter() -> Self {
        Self::new(
            "m",
            1.0,
            Dimension {
                exponents: [1, 0, 0, 0, 0, 0, 0],
            },
        )
    }

    /// Second (s).
    #[must_use]
    pub fn second() -> Self {
        Self::new(
            "s",
            1.0,
            Dimension {
                exponents: [0, 0, 1, 0, 0, 0, 0],
            },
        )
    }

    /// Kilogram (kg).
    #[must_use]
    pub fn kilogram() -> Self {
        Self::new(
            "kg",
            1.0,
            Dimension {
                exponents: [0, 1, 0, 0, 0, 0, 0],
            },
        )
    }

    /// Centimeter (cm).
    #[must_use]
    pub fn centimeter() -> Self {
        Self::new(
            "cm",
            0.01,
            Dimension {
                exponents: [1, 0, 0, 0, 0, 0, 0],
            },
        )
    }

    /// Newton (kg·m/s^2).
    #[must_use]
    pub fn newton() -> Self {
        Self::new(
            "N",
            1.0,
            Dimension {
                exponents: [1, 1, -2, 0, 0, 0, 0],
            },
        )
    }

    /// Multiplies two units.
    #[must_use]
    pub fn mul(&self, other: &Self, symbol: &'static str) -> Self {
        let mut exponents = [0; 7];
        for (i, exp) in exponents.iter_mut().enumerate() {
            *exp = self.dimension.exponents[i] + other.dimension.exponents[i];
        }

        Self::new(
            symbol,
            self.scale_to_si * other.scale_to_si,
            Dimension { exponents },
        )
    }

    /// Divides two units.
    #[must_use]
    pub fn div(&self, other: &Self, symbol: &'static str) -> Self {
        let mut exponents = [0; 7];
        for (i, exp) in exponents.iter_mut().enumerate() {
            *exp = self.dimension.exponents[i] - other.dimension.exponents[i];
        }

        Self::new(
            symbol,
            self.scale_to_si / other.scale_to_si,
            Dimension { exponents },
        )
    }

    /// Raises a unit to an integer power.
    #[must_use]
    pub fn powi(&self, power: i8, symbol: &'static str) -> Self {
        Self::new(
            symbol,
            self.scale_to_si.powi(i32::from(power)),
            self.dimension.powi(power),
        )
    }

    /// Returns true if units have compatible dimensions.
    #[must_use]
    pub fn is_compatible_with(&self, other: &Self) -> bool {
        self.dimension == other.dimension
    }

    /// Returns scale factor to convert from self to other.
    ///
    /// Returns `None` when dimensions are incompatible.
    #[must_use]
    pub fn conversion_factor_to(&self, other: &Self) -> Option<f64> {
        if !self.is_compatible_with(other) {
            return None;
        }
        Some(self.scale_to_si / other.scale_to_si)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dimension_multiplication_and_division() {
        let length = Unit::meter();
        let time = Unit::second();
        let velocity = length.div(&time, "m/s");
        let acceleration = velocity.div(&time, "m/s^2");

        assert_eq!(velocity.dimension.exponents, [1, 0, -1, 0, 0, 0, 0]);
        assert_eq!(acceleration.dimension.exponents, [1, 0, -2, 0, 0, 0, 0]);
    }

    #[test]
    fn test_newton_dimension() {
        let n = Unit::newton();
        assert_eq!(n.dimension.exponents, [1, 1, -2, 0, 0, 0, 0]);
    }

    #[test]
    fn test_unit_compatibility_and_conversion() {
        let m = Unit::meter();
        let cm = Unit::centimeter();
        assert!(cm.is_compatible_with(&m));

        let factor = cm.conversion_factor_to(&m).unwrap();
        assert!((factor - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_incompatible_units() {
        let m = Unit::meter();
        let s = Unit::second();
        assert!(!m.is_compatible_with(&s));
        assert_eq!(m.conversion_factor_to(&s), None);
    }
}
