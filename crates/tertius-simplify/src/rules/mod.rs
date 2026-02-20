//! Simplification rules organized by category.

pub mod arithmetic;
pub mod exp_log;
pub mod integration;
pub mod trig;

use crate::language::TertiusLang;
use egg::Rewrite;

/// Collects all simplification rules.
pub fn all_rules() -> Vec<Rewrite<TertiusLang, ()>> {
    let mut rules = Vec::new();
    rules.extend(arithmetic::rules());
    rules.extend(trig::rules());
    rules.extend(exp_log::rules_branch_safe());
    rules.extend(integration::rules());
    rules
}

/// Collects domain-sensitive rules that should be enabled only under
/// explicit assumptions (for example, all symbols real).
pub fn real_assumption_rules() -> Vec<Rewrite<TertiusLang, ()>> {
    exp_log::rules_real_domain()
}
