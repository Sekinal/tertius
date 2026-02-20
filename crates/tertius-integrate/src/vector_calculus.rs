//! Symbolic scaffolding for vector-calculus integrals.

use tertius_core::arena::ExprArena;
use tertius_core::expr::ExprNode;
use tertius_core::handle::ExprHandle;

/// Coordinate frame for vector-calculus integrals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoordinateFrame {
    /// `(x, y, z)`
    Cartesian,
    /// `(r, theta, z)`
    Cylindrical,
    /// `(r, theta, phi)`
    Spherical,
}

/// Line integral scaffold `∫ F · dr`.
#[derive(Clone, Debug)]
pub struct LineIntegralScaffold {
    /// Coordinate frame for interpretation.
    pub frame: CoordinateFrame,
    /// Integrand after substituting a parametrized tangent vector.
    pub integrand: ExprHandle,
    /// Parameter of integration.
    pub parameter: ExprHandle,
    /// Lower parameter bound.
    pub lower: ExprHandle,
    /// Upper parameter bound.
    pub upper: ExprHandle,
}

/// Surface flux integral scaffold `∬ F · dS`.
#[derive(Clone, Debug)]
pub struct SurfaceIntegralScaffold {
    /// Coordinate frame for interpretation.
    pub frame: CoordinateFrame,
    /// Integrand `F · (n dS)`.
    pub integrand: ExprHandle,
    /// First surface parameter.
    pub u: ExprHandle,
    /// Second surface parameter.
    pub v: ExprHandle,
    /// Lower bound for `u`.
    pub u_lower: ExprHandle,
    /// Upper bound for `u`.
    pub u_upper: ExprHandle,
    /// Lower bound for `v`.
    pub v_lower: ExprHandle,
    /// Upper bound for `v`.
    pub v_upper: ExprHandle,
}

/// Volume integral scaffold `∭ f J dV`.
#[derive(Clone, Debug)]
pub struct VolumeIntegralScaffold {
    /// Coordinate frame for interpretation.
    pub frame: CoordinateFrame,
    /// Integrand `f * jacobian`.
    pub integrand: ExprHandle,
    /// Variables of integration.
    pub vars: [ExprHandle; 3],
    /// Lower bounds.
    pub lower: [ExprHandle; 3],
    /// Upper bounds.
    pub upper: [ExprHandle; 3],
}

/// Build a line integral scaffold from field and tangent.
#[must_use]
pub fn line_integral_from_tangent(
    arena: &mut ExprArena,
    frame: CoordinateFrame,
    field: [ExprHandle; 3],
    tangent: [ExprHandle; 3],
    parameter: ExprHandle,
    lower: ExprHandle,
    upper: ExprHandle,
) -> LineIntegralScaffold {
    LineIntegralScaffold {
        frame,
        integrand: dot3(arena, field, tangent),
        parameter,
        lower,
        upper,
    }
}

/// Build a surface flux scaffold from field and oriented area vector.
#[must_use]
pub fn surface_flux_integral(
    arena: &mut ExprArena,
    frame: CoordinateFrame,
    field: [ExprHandle; 3],
    oriented_area: [ExprHandle; 3],
    u: ExprHandle,
    v: ExprHandle,
    u_lower: ExprHandle,
    u_upper: ExprHandle,
    v_lower: ExprHandle,
    v_upper: ExprHandle,
) -> SurfaceIntegralScaffold {
    SurfaceIntegralScaffold {
        frame,
        integrand: dot3(arena, field, oriented_area),
        u,
        v,
        u_lower,
        u_upper,
        v_lower,
        v_upper,
    }
}

/// Build a volume integral scaffold for scalar density with jacobian.
#[must_use]
pub fn volume_integral(
    arena: &mut ExprArena,
    frame: CoordinateFrame,
    scalar_density: ExprHandle,
    jacobian: ExprHandle,
    vars: [ExprHandle; 3],
    lower: [ExprHandle; 3],
    upper: [ExprHandle; 3],
) -> VolumeIntegralScaffold {
    let integrand = arena.mul(smallvec::smallvec![scalar_density, jacobian]);
    VolumeIntegralScaffold {
        frame,
        integrand,
        vars,
        lower,
        upper,
    }
}

fn dot3(arena: &mut ExprArena, a: [ExprHandle; 3], b: [ExprHandle; 3]) -> ExprHandle {
    let p0 = arena.mul(smallvec::smallvec![a[0], b[0]]);
    let p1 = arena.mul(smallvec::smallvec![a[1], b[1]]);
    let p2 = arena.mul(smallvec::smallvec![a[2], b[2]]);
    arena.intern(ExprNode::Add(smallvec::smallvec![p0, p1, p2]))
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use smallvec::smallvec;
    use tertius_core::expr::{functions, ExprNode};

    use super::*;

    fn eval(arena: &ExprArena, expr: ExprHandle, vars: &HashMap<&str, f64>) -> f64 {
        match arena.get(expr) {
            ExprNode::Integer(n) => *n as f64,
            ExprNode::Rational(num, den) => *num as f64 / *den as f64,
            ExprNode::Symbol(id) => {
                let name = arena.symbol_name(*id).unwrap();
                *vars.get(name).unwrap()
            }
            ExprNode::Add(args) => args.iter().map(|h| eval(arena, *h, vars)).sum(),
            ExprNode::Mul(args) => args.iter().map(|h| eval(arena, *h, vars)).product(),
            ExprNode::Pow { base, exp } => eval(arena, *base, vars).powf(eval(arena, *exp, vars)),
            ExprNode::Neg(arg) => -eval(arena, *arg, vars),
            ExprNode::Div { num, den } => eval(arena, *num, vars) / eval(arena, *den, vars),
            ExprNode::Function { id, args } => {
                let x = eval(arena, args[0], vars);
                match *id {
                    functions::SIN => x.sin(),
                    functions::COS => x.cos(),
                    functions::EXP => x.exp(),
                    _ => panic!("unsupported function"),
                }
            }
        }
    }

    #[test]
    fn test_line_integral_scaffold_uses_dot_product() {
        let mut arena = ExprArena::new();
        let t = arena.symbol("t");
        let one = arena.integer(1);
        let two = arena.integer(2);
        let zero_bound = arena.integer(0);
        let one_bound = arena.integer(1);
        let field = [t, one, two];
        let tangent = [one, two, t];
        let li = line_integral_from_tangent(
            &mut arena,
            CoordinateFrame::Cartesian,
            field,
            tangent,
            t,
            zero_bound,
            one_bound,
        );
        // dot = t*1 + 1*2 + 2*t = 3t + 2
        let vars = HashMap::from([("t", 4.0)]);
        assert!((eval(&arena, li.integrand, &vars) - 14.0).abs() < 1e-9);
    }

    #[test]
    fn test_surface_flux_scaffold_uses_dot_product() {
        let mut arena = ExprArena::new();
        let u = arena.symbol("u");
        let v = arena.symbol("v");
        let zero_bound = arena.integer(0);
        let one_bound = arena.integer(1);
        let field = [u, v, arena.integer(1)];
        let area = [arena.integer(0), arena.integer(0), arena.integer(3)];
        let s = surface_flux_integral(
            &mut arena,
            CoordinateFrame::Cartesian,
            field,
            area,
            u,
            v,
            zero_bound,
            one_bound,
            zero_bound,
            one_bound,
        );
        let vars = HashMap::from([("u", 2.0), ("v", -1.0)]);
        assert!((eval(&arena, s.integrand, &vars) - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_volume_integral_scaffold_applies_jacobian() {
        let mut arena = ExprArena::new();
        let r = arena.symbol("r");
        let theta = arena.symbol("theta");
        let z = arena.symbol("z");
        let density = arena.add(smallvec![r, z]);
        let jac = r;
        let lower = [arena.integer(0), arena.integer(0), arena.integer(0)];
        let upper = [arena.integer(1), arena.integer(1), arena.integer(1)];
        let v = volume_integral(
            &mut arena,
            CoordinateFrame::Cylindrical,
            density,
            jac,
            [r, theta, z],
            lower,
            upper,
        );
        // (r+z) * r
        let vars = HashMap::from([("r", 2.0), ("theta", 0.5), ("z", 3.0)]);
        assert!((eval(&arena, v.integrand, &vars) - 10.0).abs() < 1e-9);
    }
}
