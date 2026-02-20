//! Physics-focused facade APIs.

pub mod dashboard;
pub mod latex;
pub mod mode;
pub mod trace;

pub use dashboard::{CorrectnessReport, PerformanceReport, PhysicsDashboard};
pub use latex::{expr_to_latex, format_quantity_latex, parse_linear_latex, parse_quantity_latex};
pub use mode::{FirstOrderLinearRequest, PhysicsMode, PhysicsSolveResult};
pub use trace::{DerivationStep, DerivationTrace};
