//! Coordinate-aware vector calculus operators.

use smallvec::smallvec;
use tertius_core::arena::ExprArena;
use tertius_core::expr::{functions, ExprNode};
use tertius_core::handle::ExprHandle;

use crate::partial::partial;

/// Gradient in Cartesian coordinates.
pub fn gradient_cartesian(
    arena: &mut ExprArena,
    scalar: ExprHandle,
    vars_xyz: [ExprHandle; 3],
) -> [ExprHandle; 3] {
    [
        partial(arena, scalar, vars_xyz[0]),
        partial(arena, scalar, vars_xyz[1]),
        partial(arena, scalar, vars_xyz[2]),
    ]
}

/// Divergence in Cartesian coordinates.
pub fn divergence_cartesian(
    arena: &mut ExprArena,
    vector_xyz: [ExprHandle; 3],
    vars_xyz: [ExprHandle; 3],
) -> ExprHandle {
    let dfx_dx = partial(arena, vector_xyz[0], vars_xyz[0]);
    let dfy_dy = partial(arena, vector_xyz[1], vars_xyz[1]);
    let dfz_dz = partial(arena, vector_xyz[2], vars_xyz[2]);
    arena.add(smallvec![dfx_dx, dfy_dy, dfz_dz])
}

/// Curl in Cartesian coordinates.
pub fn curl_cartesian(
    arena: &mut ExprArena,
    vector_xyz: [ExprHandle; 3],
    vars_xyz: [ExprHandle; 3],
) -> [ExprHandle; 3] {
    let [fx, fy, fz] = vector_xyz;
    let [x, y, z] = vars_xyz;

    let dfz_dy = partial(arena, fz, y);
    let dfy_dz = partial(arena, fy, z);
    let c1 = sub(arena, dfz_dy, dfy_dz);

    let dfx_dz = partial(arena, fx, z);
    let dfz_dx = partial(arena, fz, x);
    let c2 = sub(arena, dfx_dz, dfz_dx);

    let dfy_dx = partial(arena, fy, x);
    let dfx_dy = partial(arena, fx, y);
    let c3 = sub(arena, dfy_dx, dfx_dy);

    [c1, c2, c3]
}

/// Gradient in cylindrical coordinates (r, theta, z).
pub fn gradient_cylindrical(
    arena: &mut ExprArena,
    scalar: ExprHandle,
    r: ExprHandle,
    theta: ExprHandle,
    z: ExprHandle,
) -> [ExprHandle; 3] {
    let df_dr = partial(arena, scalar, r);
    let df_dtheta = partial(arena, scalar, theta);
    let df_dz = partial(arena, scalar, z);
    [df_dr, div(arena, df_dtheta, r), df_dz]
}

/// Gradient in spherical coordinates (r, theta, phi).
pub fn gradient_spherical(
    arena: &mut ExprArena,
    scalar: ExprHandle,
    r: ExprHandle,
    theta: ExprHandle,
    phi: ExprHandle,
) -> [ExprHandle; 3] {
    let df_dr = partial(arena, scalar, r);
    let df_dtheta = partial(arena, scalar, theta);
    let df_dphi = partial(arena, scalar, phi);

    let sin_theta = arena.intern(ExprNode::Function {
        id: functions::SIN,
        args: smallvec![theta],
    });
    let r_sin_theta = arena.mul(smallvec![r, sin_theta]);
    [
        df_dr,
        div(arena, df_dtheta, r),
        div(arena, df_dphi, r_sin_theta),
    ]
}

fn sub(arena: &mut ExprArena, a: ExprHandle, b: ExprHandle) -> ExprHandle {
    let neg_b = arena.neg(b);
    arena.add(smallvec![a, neg_b])
}

fn div(arena: &mut ExprArena, num: ExprHandle, den: ExprHandle) -> ExprHandle {
    arena.intern(ExprNode::Div { num, den })
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
                    functions::TAN => x.tan(),
                    functions::EXP => x.exp(),
                    functions::LN => x.ln(),
                    functions::SQRT => x.sqrt(),
                    functions::ABS => x.abs(),
                    _ => panic!("unsupported function"),
                }
            }
        }
    }

    #[test]
    fn test_gradient_cartesian_quadratic() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let y = arena.symbol("y");
        let z = arena.symbol("z");
        let two = arena.integer(2);
        let x2 = arena.pow(x, two);
        let y2 = arena.pow(y, two);
        let z2 = arena.pow(z, two);
        let f = arena.add(smallvec![x2, y2, z2]);

        let g = gradient_cartesian(&mut arena, f, [x, y, z]);
        let vars = HashMap::from([("x", 1.0), ("y", 2.0), ("z", 3.0)]);
        assert!((eval(&arena, g[0], &vars) - 2.0).abs() < 1e-9);
        assert!((eval(&arena, g[1], &vars) - 4.0).abs() < 1e-9);
        assert!((eval(&arena, g[2], &vars) - 6.0).abs() < 1e-9);
    }

    #[test]
    fn test_divergence_cartesian_identity_field() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let y = arena.symbol("y");
        let z = arena.symbol("z");

        let div = divergence_cartesian(&mut arena, [x, y, z], [x, y, z]);
        let vars = HashMap::from([("x", 4.0), ("y", -1.0), ("z", 2.0)]);
        assert!((eval(&arena, div, &vars) - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_curl_cartesian_rotation_field() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let y = arena.symbol("y");
        let z = arena.symbol("z");
        let zero = arena.integer(0);
        let neg_y = arena.neg(y);
        let field = [neg_y, x, zero];

        let c = curl_cartesian(&mut arena, field, [x, y, z]);
        let vars = HashMap::from([("x", 1.0), ("y", 2.0), ("z", 3.0)]);
        assert!(eval(&arena, c[0], &vars).abs() < 1e-9);
        assert!(eval(&arena, c[1], &vars).abs() < 1e-9);
        assert!((eval(&arena, c[2], &vars) - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_gradient_cylindrical_radial() {
        let mut arena = ExprArena::new();
        let r = arena.symbol("r");
        let theta = arena.symbol("theta");
        let z = arena.symbol("z");
        let two = arena.integer(2);
        let r2 = arena.pow(r, two);
        let z2 = arena.pow(z, two);
        let f = arena.add(smallvec![r2, z2]);

        let g = gradient_cylindrical(&mut arena, f, r, theta, z);
        let vars = HashMap::from([("r", 3.0), ("theta", 0.7), ("z", 4.0)]);
        assert!((eval(&arena, g[0], &vars) - 6.0).abs() < 1e-9);
        assert!(eval(&arena, g[1], &vars).abs() < 1e-9);
        assert!((eval(&arena, g[2], &vars) - 8.0).abs() < 1e-9);
    }

    #[test]
    fn test_gradient_spherical_radius() {
        let mut arena = ExprArena::new();
        let r = arena.symbol("r");
        let theta = arena.symbol("theta");
        let phi = arena.symbol("phi");

        let g = gradient_spherical(&mut arena, r, r, theta, phi);
        let vars = HashMap::from([("r", 2.0), ("theta", 1.0), ("phi", 0.5)]);
        assert!((eval(&arena, g[0], &vars) - 1.0).abs() < 1e-9);
        assert!(eval(&arena, g[1], &vars).abs() < 1e-9);
        assert!(eval(&arena, g[2], &vars).abs() < 1e-9);
    }
}
