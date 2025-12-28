//! Polynomial system solving for Tertius CAS.
//!
//! This crate provides algorithms for solving systems of polynomial equations:
//!
//! - **FGLM**: Basis conversion from grevlex to lex ordering
//! - **Triangular decomposition**: Extract solutions from lex Gröbner bases
//! - **RUR**: Rational Univariate Representation for algebraic solutions
//!
//! # Example
//!
//! ```ignore
//! use tertius_solve::solve_system;
//! use tertius_rings::finite_field::FiniteField;
//!
//! type GF101 = FiniteField<101>;
//!
//! // System: x + y = 2, x - y = 0 over GF(101)
//! // Solution: x = 1, y = 1
//! let system = vec![
//!     // x + y - 2
//!     vec![(GF101::new(1), mono(&[1, 0])), (GF101::new(1), mono(&[0, 1])), (GF101::new(99), mono(&[0, 0]))],
//!     // x - y
//!     vec![(GF101::new(1), mono(&[1, 0])), (GF101::new(100), mono(&[0, 1]))],
//! ];
//!
//! let solutions = solve_system(system);
//! ```

pub mod fglm;
pub mod quotient;
pub mod triangular;

pub use fglm::{fglm_convert, FGLMResult};
pub use triangular::{solve_triangular, TriangularSystem};
