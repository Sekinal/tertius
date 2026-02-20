//! Simplification rules organized by category.

pub mod arithmetic;
pub mod trig;
pub mod exp_log;
pub mod integration;

use egg::Rewrite;
use crate::language::TertiusLang;

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
