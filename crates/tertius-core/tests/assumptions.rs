use tertius_core::assumptions::{Assumption, AssumptionSet};

#[test]
fn test_assumption_set_basic_membership() {
    let mut assumptions = AssumptionSet::new();
    assumptions.assume("x", Assumption::Real);
    assumptions.assume("x", Assumption::Positive);

    assert!(assumptions.has("x", Assumption::Real));
    assert!(assumptions.has("x", Assumption::Positive));
    assert!(!assumptions.has("x", Assumption::NonZero));
    assert!(!assumptions.has("y", Assumption::Real));
}

#[test]
fn test_assumption_set_real_coverage() {
    let mut assumptions = AssumptionSet::new();
    assumptions.assume("x", Assumption::Real);
    assumptions.assume("y", Assumption::Real);

    assert!(assumptions.all_real(["x", "y"]));
    assert!(!assumptions.all_real(["x", "z"]));
}
