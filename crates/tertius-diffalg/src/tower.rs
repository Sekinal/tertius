//! Transcendental tower for the Risch algorithm.
//!
//! A transcendental tower is a sequence of field extensions:
//!
//! K = K₀ ⊂ K₁ = K₀(θ₁) ⊂ K₂ = K₁(θ₂) ⊂ ...
//!
//! where each θᵢ is either:
//! - Logarithmic: θᵢ = log(uᵢ) for some uᵢ ∈ Kᵢ₋₁
//! - Exponential: θᵢ = exp(uᵢ) for some uᵢ ∈ Kᵢ₋₁
//!
//! The tower allows us to represent expressions involving nested
//! logarithms and exponentials in a structured way.

use std::fmt;

/// The type of a transcendental extension.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TranscendentalType {
    /// θ = log(u) where u is represented at a previous level
    Logarithmic {
        /// The level where the argument u is defined
        argument_level: usize,
        /// String representation of the argument
        argument_repr: String,
    },
    /// θ = exp(u) where u is represented at a previous level
    Exponential {
        /// The level where the exponent u is defined
        exponent_level: usize,
        /// String representation of the exponent
        exponent_repr: String,
    },
}

impl TranscendentalType {
    /// Returns true if this is a logarithmic extension.
    pub fn is_logarithmic(&self) -> bool {
        matches!(self, TranscendentalType::Logarithmic { .. })
    }

    /// Returns true if this is an exponential extension.
    pub fn is_exponential(&self) -> bool {
        matches!(self, TranscendentalType::Exponential { .. })
    }
}

impl fmt::Display for TranscendentalType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TranscendentalType::Logarithmic { argument_repr, .. } => {
                write!(f, "log({})", argument_repr)
            }
            TranscendentalType::Exponential { exponent_repr, .. } => {
                write!(f, "exp({})", exponent_repr)
            }
        }
    }
}

/// A single level in the transcendental tower.
#[derive(Clone, Debug)]
pub struct TowerLevel {
    /// The name of the transcendental at this level (e.g., "θ₁")
    pub name: String,
    /// The type of extension
    pub extension_type: TranscendentalType,
    /// The level index (0-based)
    pub level: usize,
}

impl TowerLevel {
    /// Creates a new tower level.
    pub fn new(name: String, extension_type: TranscendentalType, level: usize) -> Self {
        Self {
            name,
            extension_type,
            level,
        }
    }

    /// Creates a logarithmic level: θ = log(argument)
    pub fn logarithmic(name: &str, level: usize, arg_level: usize, arg_repr: &str) -> Self {
        Self {
            name: name.to_string(),
            extension_type: TranscendentalType::Logarithmic {
                argument_level: arg_level,
                argument_repr: arg_repr.to_string(),
            },
            level,
        }
    }

    /// Creates an exponential level: θ = exp(exponent)
    pub fn exponential(name: &str, level: usize, exp_level: usize, exp_repr: &str) -> Self {
        Self {
            name: name.to_string(),
            extension_type: TranscendentalType::Exponential {
                exponent_level: exp_level,
                exponent_repr: exp_repr.to_string(),
            },
            level,
        }
    }
}

impl fmt::Display for TowerLevel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} = {}", self.name, self.extension_type)
    }
}

/// A complete transcendental tower.
///
/// Represents the field extension K(θ₁, θ₂, ..., θₙ) where each θᵢ
/// is a transcendental (logarithm or exponential) over the previous level.
#[derive(Clone, Debug)]
pub struct TranscendentalTower {
    /// The base variable name (usually "x")
    pub base_variable: String,
    /// The levels of the tower, from bottom (θ₁) to top (θₙ)
    pub levels: Vec<TowerLevel>,
}

impl TranscendentalTower {
    /// Creates a new tower with the given base variable.
    pub fn new(base_variable: &str) -> Self {
        Self {
            base_variable: base_variable.to_string(),
            levels: Vec::new(),
        }
    }

    /// Returns the height of the tower (number of transcendental extensions).
    pub fn height(&self) -> usize {
        self.levels.len()
    }

    /// Returns true if the tower has no extensions (just the base field).
    pub fn is_base_field(&self) -> bool {
        self.levels.is_empty()
    }

    /// Adds a logarithmic extension: θ = log(argument).
    pub fn add_logarithmic(&mut self, arg_level: usize, arg_repr: &str) -> &TowerLevel {
        let level = self.levels.len();
        let name = format!("θ{}", level + 1);
        self.levels
            .push(TowerLevel::logarithmic(&name, level, arg_level, arg_repr));
        self.levels.last().unwrap()
    }

    /// Adds an exponential extension: θ = exp(exponent).
    pub fn add_exponential(&mut self, exp_level: usize, exp_repr: &str) -> &TowerLevel {
        let level = self.levels.len();
        let name = format!("θ{}", level + 1);
        self.levels
            .push(TowerLevel::exponential(&name, level, exp_level, exp_repr));
        self.levels.last().unwrap()
    }

    /// Returns the level at the given index.
    pub fn level(&self, index: usize) -> Option<&TowerLevel> {
        self.levels.get(index)
    }

    /// Returns the top level of the tower.
    pub fn top(&self) -> Option<&TowerLevel> {
        self.levels.last()
    }

    /// Checks if all extensions are logarithmic (purely logarithmic tower).
    pub fn is_purely_logarithmic(&self) -> bool {
        self.levels.iter().all(|l| l.extension_type.is_logarithmic())
    }

    /// Checks if all extensions are exponential (purely exponential tower).
    pub fn is_purely_exponential(&self) -> bool {
        self.levels
            .iter()
            .all(|l| l.extension_type.is_exponential())
    }

    /// Returns a string representation of the tower.
    pub fn describe(&self) -> String {
        let mut desc = format!("K = Q({})", self.base_variable);

        for (i, level) in self.levels.iter().enumerate() {
            desc.push_str(&format!("\nK{} = K{}({}) where {}", i + 1, i, level.name, level));
        }

        desc
    }
}

impl fmt::Display for TranscendentalTower {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Q({})", self.base_variable)?;
        for level in &self.levels {
            write!(f, "({})", level.name)?;
        }
        Ok(())
    }
}

/// Common transcendental tower patterns.
impl TranscendentalTower {
    /// Creates a tower for log(x): Q(x)(θ) where θ = log(x)
    pub fn single_log(var: &str) -> Self {
        let mut tower = Self::new(var);
        tower.add_logarithmic(0, var);
        tower
    }

    /// Creates a tower for exp(x): Q(x)(θ) where θ = exp(x)
    pub fn single_exp(var: &str) -> Self {
        let mut tower = Self::new(var);
        tower.add_exponential(0, var);
        tower
    }

    /// Creates a tower for log(log(x)): Q(x)(θ₁)(θ₂) where θ₁ = log(x), θ₂ = log(θ₁)
    pub fn double_log(var: &str) -> Self {
        let mut tower = Self::new(var);
        tower.add_logarithmic(0, var);
        tower.add_logarithmic(1, "θ1");
        tower
    }

    /// Creates a tower for exp(exp(x)): Q(x)(θ₁)(θ₂) where θ₁ = exp(x), θ₂ = exp(θ₁)
    pub fn double_exp(var: &str) -> Self {
        let mut tower = Self::new(var);
        tower.add_exponential(0, var);
        tower.add_exponential(1, "θ1");
        tower
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_tower() {
        let tower = TranscendentalTower::new("x");

        assert_eq!(tower.height(), 0);
        assert!(tower.is_base_field());
        assert!(tower.is_purely_logarithmic()); // vacuously true
        assert!(tower.is_purely_exponential()); // vacuously true
    }

    #[test]
    fn test_single_log_tower() {
        let tower = TranscendentalTower::single_log("x");

        assert_eq!(tower.height(), 1);
        assert!(!tower.is_base_field());
        assert!(tower.is_purely_logarithmic());
        assert!(!tower.is_purely_exponential());

        let top = tower.top().unwrap();
        assert!(top.extension_type.is_logarithmic());
    }

    #[test]
    fn test_single_exp_tower() {
        let tower = TranscendentalTower::single_exp("x");

        assert_eq!(tower.height(), 1);
        assert!(tower.is_purely_exponential());
        assert!(!tower.is_purely_logarithmic());
    }

    #[test]
    fn test_mixed_tower() {
        let mut tower = TranscendentalTower::new("x");
        tower.add_logarithmic(0, "x");
        tower.add_exponential(1, "θ1");

        assert_eq!(tower.height(), 2);
        assert!(!tower.is_purely_logarithmic());
        assert!(!tower.is_purely_exponential());
    }

    #[test]
    fn test_tower_display() {
        let tower = TranscendentalTower::single_log("x");
        let s = format!("{}", tower);
        assert!(s.contains("θ1"));
    }

    #[test]
    fn test_tower_describe() {
        let tower = TranscendentalTower::single_log("x");
        let desc = tower.describe();

        assert!(desc.contains("Q(x)"));
        assert!(desc.contains("log"));
    }
}
