//! Differentiation rules documentation.
//!
//! This module documents the differentiation rules implemented by this crate.
//!
//! # Basic Rules
//!
//! | Rule | Formula |
//! |------|---------|
//! | Constant | d/dx(c) = 0 |
//! | Variable | d/dx(x) = 1 |
//! | Sum | d/dx(f + g) = f' + g' |
//! | Product | d/dx(f · g) = f'g + fg' |
//! | Quotient | d/dx(f/g) = (f'g - fg')/g² |
//! | Power | d/dx(xⁿ) = n·xⁿ⁻¹ |
//! | Chain | d/dx(f(g)) = f'(g) · g' |
//!
//! # Elementary Functions
//!
//! | Function | Derivative |
//! |----------|------------|
//! | sin(x) | cos(x) |
//! | cos(x) | -sin(x) |
//! | tan(x) | sec²(x) = 1/cos²(x) |
//! | exp(x) | exp(x) |
//! | ln(x) | 1/x |
//! | log₁₀(x) | 1/(x·ln(10)) |
//! | √x | 1/(2√x) |
//! | \|x\| | x/\|x\| (sign function) |
//!
//! # Generalized Power Rule
//!
//! For f^g where both f and g depend on x:
//!
//! d/dx(f^g) = f^g · (g' · ln(f) + g · f'/f)
//!
//! This formula is derived using logarithmic differentiation.
