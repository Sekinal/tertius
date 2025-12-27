//! Simplification rules organized by category.

pub mod arithmetic;
pub mod trig;
pub mod exp_log;

use egg::Rewrite;
use crate::language::TertiusLang;

/// Collects all simplification rules.
pub fn all_rules() -> Vec<Rewrite<TertiusLang, ()>> {
    let mut rules = Vec::new();
    rules.extend(arithmetic::rules());
    rules.extend(trig::rules());
    rules.extend(exp_log::rules());
    rules
}
