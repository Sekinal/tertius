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

/// Divergence in cylindrical coordinates (r, theta, z).
pub fn divergence_cylindrical(
    arena: &mut ExprArena,
    vector_rtz: [ExprHandle; 3],
    r: ExprHandle,
    theta: ExprHandle,
    z: ExprHandle,
) -> ExprHandle {
    let [a_r, a_theta, a_z] = vector_rtz;
    let r_a_r = mul(arena, r, a_r);
    let d_r_ar_dr = partial(arena, r_a_r, r);
    let d_atheta_dtheta = partial(arena, a_theta, theta);
    let d_az_dz = partial(arena, a_z, z);

    let term_r = div(arena, d_r_ar_dr, r);
    let term_theta = div(arena, d_atheta_dtheta, r);
    arena.add(smallvec![term_r, term_theta, d_az_dz])
}

/// Curl in cylindrical coordinates (r, theta, z).
pub fn curl_cylindrical(
    arena: &mut ExprArena,
    vector_rtz: [ExprHandle; 3],
    r: ExprHandle,
    theta: ExprHandle,
    z: ExprHandle,
) -> [ExprHandle; 3] {
    let [a_r, a_theta, a_z] = vector_rtz;

    let d_az_dtheta = partial(arena, a_z, theta);
    let d_atheta_dz = partial(arena, a_theta, z);
    let scaled_d_az_dtheta = div(arena, d_az_dtheta, r);
    let c_r = sub(arena, scaled_d_az_dtheta, d_atheta_dz);

    let d_ar_dz = partial(arena, a_r, z);
    let d_az_dr = partial(arena, a_z, r);
    let c_theta = sub(arena, d_ar_dz, d_az_dr);

    let r_a_theta = mul(arena, r, a_theta);
    let d_r_atheta_dr = partial(arena, r_a_theta, r);
    let d_ar_dtheta = partial(arena, a_r, theta);
    let z_num = sub(arena, d_r_atheta_dr, d_ar_dtheta);
    let c_z = div(arena, z_num, r);

    [c_r, c_theta, c_z]
}

/// Divergence in spherical coordinates (r, theta, phi).
pub fn divergence_spherical(
    arena: &mut ExprArena,
    vector_rtp: [ExprHandle; 3],
    r: ExprHandle,
    theta: ExprHandle,
    phi: ExprHandle,
) -> ExprHandle {
    let [a_r, a_theta, a_phi] = vector_rtp;
    let two = arena.integer(2);
    let sin_theta = arena.intern(ExprNode::Function {
        id: functions::SIN,
        args: smallvec![theta],
    });
    let r2 = arena.pow(r, two);
    let r2_ar = mul(arena, r2, a_r);
    let d_r2ar_dr = partial(arena, r2_ar, r);
    let term_r = div(arena, d_r2ar_dr, r2);

    let sin_a_theta = mul(arena, sin_theta, a_theta);
    let d_sin_atheta_dtheta = partial(arena, sin_a_theta, theta);
    let r_sin_theta = mul(arena, r, sin_theta);
    let term_theta = div(arena, d_sin_atheta_dtheta, r_sin_theta);

    let d_aphi_dphi = partial(arena, a_phi, phi);
    let term_phi = div(arena, d_aphi_dphi, r_sin_theta);
    arena.add(smallvec![term_r, term_theta, term_phi])
}

/// Curl in spherical coordinates (r, theta, phi).
pub fn curl_spherical(
    arena: &mut ExprArena,
    vector_rtp: [ExprHandle; 3],
    r: ExprHandle,
    theta: ExprHandle,
    phi: ExprHandle,
) -> [ExprHandle; 3] {
    let [a_r, a_theta, a_phi] = vector_rtp;
    let sin_theta = arena.intern(ExprNode::Function {
        id: functions::SIN,
        args: smallvec![theta],
    });
    let r_sin_theta = mul(arena, r, sin_theta);

    let sin_a_phi = mul(arena, sin_theta, a_phi);
    let d_sin_aphi_dtheta = partial(arena, sin_a_phi, theta);
    let d_atheta_dphi = partial(arena, a_theta, phi);
    let r_num = sub(arena, d_sin_aphi_dtheta, d_atheta_dphi);
    let c_r = div(arena, r_num, r_sin_theta);

    let d_ar_dphi = partial(arena, a_r, phi);
    let left = div(arena, d_ar_dphi, sin_theta);
    let r_a_phi = mul(arena, r, a_phi);
    let d_r_aphi_dr = partial(arena, r_a_phi, r);
    let theta_num = sub(arena, left, d_r_aphi_dr);
    let c_theta = div(arena, theta_num, r);

    let r_a_theta = mul(arena, r, a_theta);
    let d_r_atheta_dr = partial(arena, r_a_theta, r);
    let d_ar_dtheta = partial(arena, a_r, theta);
    let phi_num = sub(arena, d_r_atheta_dr, d_ar_dtheta);
    let c_phi = div(arena, phi_num, r);

    [c_r, c_theta, c_phi]
}

/// Scalar Laplacian in Cartesian coordinates.
pub fn laplacian_cartesian(
    arena: &mut ExprArena,
    scalar: ExprHandle,
    vars_xyz: [ExprHandle; 3],
) -> ExprHandle {
    let dx = partial(arena, scalar, vars_xyz[0]);
    let dy = partial(arena, scalar, vars_xyz[1]);
    let dz = partial(arena, scalar, vars_xyz[2]);
    let d2x = partial(arena, dx, vars_xyz[0]);
    let d2y = partial(arena, dy, vars_xyz[1]);
    let d2z = partial(arena, dz, vars_xyz[2]);
    arena.add(smallvec![d2x, d2y, d2z])
}

/// Scalar Laplacian in cylindrical coordinates (r, theta, z).
pub fn laplacian_cylindrical(
    arena: &mut ExprArena,
    scalar: ExprHandle,
    r: ExprHandle,
    theta: ExprHandle,
    z: ExprHandle,
) -> ExprHandle {
    let df_dr = partial(arena, scalar, r);
    let r_df_dr = mul(arena, r, df_dr);
    let d_r_df_dr = partial(arena, r_df_dr, r);
    let radial = div(arena, d_r_df_dr, r);

    let dtheta = partial(arena, scalar, theta);
    let d2theta = partial(arena, dtheta, theta);
    let two = arena.integer(2);
    let r2 = arena.pow(r, two);
    let azimuth = div(arena, d2theta, r2);

    let dz = partial(arena, scalar, z);
    let d2z = partial(arena, dz, z);
    arena.add(smallvec![radial, azimuth, d2z])
}

/// Scalar Laplacian in spherical coordinates (r, theta, phi).
pub fn laplacian_spherical(
    arena: &mut ExprArena,
    scalar: ExprHandle,
    r: ExprHandle,
    theta: ExprHandle,
    phi: ExprHandle,
) -> ExprHandle {
    let two = arena.integer(2);
    let sin_theta = arena.intern(ExprNode::Function {
        id: functions::SIN,
        args: smallvec![theta],
    });
    let r2 = arena.pow(r, two);

    let df_dr = partial(arena, scalar, r);
    let r2_df_dr = mul(arena, r2, df_dr);
    let d_r2_df_dr = partial(arena, r2_df_dr, r);
    let radial = div(arena, d_r2_df_dr, r2);

    let df_dtheta = partial(arena, scalar, theta);
    let sin_df_dtheta = mul(arena, sin_theta, df_dtheta);
    let theta_term_num = partial(arena, sin_df_dtheta, theta);
    let r2_sin_theta = mul(arena, r2, sin_theta);
    let polar = div(arena, theta_term_num, r2_sin_theta);

    let dphi = partial(arena, scalar, phi);
    let d2phi = partial(arena, dphi, phi);
    let sin2 = arena.pow(sin_theta, two);
    let r2_sin2 = mul(arena, r2, sin2);
    let azimuth = div(arena, d2phi, r2_sin2);

    arena.add(smallvec![radial, polar, azimuth])
}

fn sub(arena: &mut ExprArena, a: ExprHandle, b: ExprHandle) -> ExprHandle {
    let neg_b = arena.neg(b);
    arena.add(smallvec![a, neg_b])
}

fn div(arena: &mut ExprArena, num: ExprHandle, den: ExprHandle) -> ExprHandle {
    arena.intern(ExprNode::Div { num, den })
}

fn mul(arena: &mut ExprArena, a: ExprHandle, b: ExprHandle) -> ExprHandle {
    arena.mul(smallvec![a, b])
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

    #[test]
    fn test_divergence_cylindrical_radial_identity() {
        let mut arena = ExprArena::new();
        let r = arena.symbol("r");
        let theta = arena.symbol("theta");
        let z = arena.symbol("z");
        let zero = arena.integer(0);

        let div = divergence_cylindrical(&mut arena, [r, zero, zero], r, theta, z);
        let vars = HashMap::from([("r", 2.0), ("theta", 0.4), ("z", 1.0)]);
        assert!((eval(&arena, div, &vars) - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_curl_cylindrical_azimuthal_field() {
        let mut arena = ExprArena::new();
        let r = arena.symbol("r");
        let theta = arena.symbol("theta");
        let z = arena.symbol("z");
        let zero = arena.integer(0);
        let field = [zero, r, zero];

        let curl = curl_cylindrical(&mut arena, field, r, theta, z);
        let vars = HashMap::from([("r", 3.0), ("theta", 1.1), ("z", -2.0)]);
        assert!(eval(&arena, curl[0], &vars).abs() < 1e-9);
        assert!(eval(&arena, curl[1], &vars).abs() < 1e-9);
        assert!((eval(&arena, curl[2], &vars) - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_divergence_spherical_radial_identity() {
        let mut arena = ExprArena::new();
        let r = arena.symbol("r");
        let theta = arena.symbol("theta");
        let phi = arena.symbol("phi");
        let zero = arena.integer(0);

        let div = divergence_spherical(&mut arena, [r, zero, zero], r, theta, phi);
        let vars = HashMap::from([("r", 2.0), ("theta", 0.8), ("phi", 0.2)]);
        assert!((eval(&arena, div, &vars) - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_curl_spherical_azimuthal_field() {
        let mut arena = ExprArena::new();
        let r = arena.symbol("r");
        let theta = arena.symbol("theta");
        let phi = arena.symbol("phi");
        let zero = arena.integer(0);
        let sin_theta = arena.intern(ExprNode::Function {
            id: functions::SIN,
            args: smallvec![theta],
        });
        let a_phi = arena.mul(smallvec![r, sin_theta]);
        let field = [zero, zero, a_phi];

        let curl = curl_spherical(&mut arena, field, r, theta, phi);
        let vars = HashMap::from([("r", 3.0), ("theta", 1.0), ("phi", 0.5)]);
        assert!((eval(&arena, curl[0], &vars) - 2.0 * 1.0f64.cos()).abs() < 1e-8);
        assert!((eval(&arena, curl[1], &vars) + 2.0 * 1.0f64.sin()).abs() < 1e-8);
        assert!(eval(&arena, curl[2], &vars).abs() < 1e-8);
    }

    #[test]
    fn test_laplacian_cartesian_quadratic() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let y = arena.symbol("y");
        let z = arena.symbol("z");
        let two = arena.integer(2);
        let x2 = arena.pow(x, two);
        let y2 = arena.pow(y, two);
        let z2 = arena.pow(z, two);
        let f = arena.add(smallvec![x2, y2, z2]);

        let lap = laplacian_cartesian(&mut arena, f, [x, y, z]);
        let vars = HashMap::from([("x", 1.0), ("y", 2.0), ("z", 3.0)]);
        assert!((eval(&arena, lap, &vars) - 6.0).abs() < 1e-9);
    }

    #[test]
    fn test_laplacian_cylindrical_radial_plus_axial() {
        let mut arena = ExprArena::new();
        let r = arena.symbol("r");
        let theta = arena.symbol("theta");
        let z = arena.symbol("z");
        let two = arena.integer(2);
        let r2 = arena.pow(r, two);
        let z2 = arena.pow(z, two);
        let f = arena.add(smallvec![r2, z2]);

        let lap = laplacian_cylindrical(&mut arena, f, r, theta, z);
        let vars = HashMap::from([("r", 4.0), ("theta", 0.2), ("z", 1.5)]);
        assert!((eval(&arena, lap, &vars) - 6.0).abs() < 1e-9);
    }

    #[test]
    fn test_laplacian_spherical_radius() {
        let mut arena = ExprArena::new();
        let r = arena.symbol("r");
        let theta = arena.symbol("theta");
        let phi = arena.symbol("phi");

        let lap = laplacian_spherical(&mut arena, r, r, theta, phi);
        let vars = HashMap::from([("r", 2.0), ("theta", 1.2), ("phi", 0.3)]);
        assert!((eval(&arena, lap, &vars) - 1.0).abs() < 1e-9);
    }
}
