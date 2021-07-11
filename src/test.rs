

use super::*;
use algebraic_robots::*;
use nalgebra::{ Unit, Vector3, Matrix4, DMatrix};
use screw_chains::{UniversalRobotsUR5, Revolute3Prismatic, ScrewChainLike};
use std::f32::consts::PI;

#[test]
fn test_to_algebra() {
    let twist = Twist::from_axis_angle_and_velocities(
        &AxisAngleRotation {
            axis: Unit::new_normalize(Vector3::new(1.0, 1.0, 1.0)).into_inner(),
            angle: 1.0
        },
        &Vector3::new(20.0, 30.0, 40.0)
    );
    let angle = 1.0 / 3.0f32.sqrt();
    let algebra = twist.to_algebra();
    let expected_algebra =  Matrix4::<f32>::from_row_slice(&[
              0.0, -angle,  angle, 20.0,
            angle,    0.0, -angle, 30.0,
           -angle,  angle,    0.0, 40.0,
              0.0,    0.0,    0.0,  0.0,
    ]);
    assert_eq!(algebra, expected_algebra);
}

#[test]
fn test_to_axis_angle_rotation() {
    let expected_angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(1.0, 4.0, -3.0)).into_inner(),
        angle: 0.5
    };
    let twist = Twist::from_axis_angle_and_velocities(
        &expected_angle_axis_rotation,
        &Vector3::new(20.0, 30.0, 40.0)
    );
    let axis_angle_rotation : AxisAngleRotation = twist.to_axis_angle_rotation();
    assert_eq!(axis_angle_rotation, expected_angle_axis_rotation);
}

#[test]
fn test_to_twist() {
    let angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(1.0, 4.0, -3.0)).into_inner(),
        angle: 0.5
    };
    let expected_twist = Twist::from_axis_angle_and_velocities(
        &angle_axis_rotation,
        &Vector3::new(20.0, 30.0, 40.0)
    );
    let algebra = expected_twist.to_algebra();
    let twist = algebra.to_twist();
    assert_eq!(twist, expected_twist);
}

#[test]
fn test_exp_only_translation() {
    let angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(1.0, 4.0, -3.0)).into_inner(),
        angle: 0.0
    };
    let twist = Twist::from_axis_angle_and_velocities(
        &angle_axis_rotation,
        &Vector3::new(20.0, -30.0, 40.0)
    );
    let algebra = twist.to_algebra();
    let transformation = algebra.exponential();
    let expected_transformation = Matrix4::<f32>::from_row_slice(&[
        1.0, 0.0, 0.0,  20.0,
        0.0, 1.0, 0.0, -30.0,
        0.0, 0.0, 1.0,  40.0,
        0.0, 0.0, 0.0,   1.0,
    ]);
    assert_eq!(transformation, expected_transformation);
}

#[test]
fn test_exp_only_rotation_x() {
    let angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(1.0, 0.0, 0.0)).into_inner(),
        angle: 1.0
    };
    let twist = Twist::from_axis_angle_and_velocities(
        &angle_axis_rotation,
        &Vector3::new(0.0, 0.0, 0.0)
    );
    let algebra = twist.to_algebra();
    let transformation = algebra.exponential();
    let expected_transformation = Matrix4::<f32>::from_row_slice(&[
        1.0, 0.0, 0.0, 0.0,
        0.0, angle_axis_rotation.angle.cos(), -angle_axis_rotation.angle.sin(), 0.0,
        0.0, angle_axis_rotation.angle.sin(), angle_axis_rotation.angle.cos(), 0.0,
        0.0, 0.0, 0.0, 1.0                            ,
    ]);
    let errors = &( transformation - expected_transformation );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_exp_only_rotation_y() {
    let angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(0.0, 1.0, 0.0)).into_inner(),
        angle: 2.0
    };
    let twist = Twist::from_axis_angle_and_velocities(
        &angle_axis_rotation,
        &Vector3::new(0.0, 0.0, 0.0)
    );
    let algebra = twist.to_algebra();
    let transformation = algebra.exponential();
    let expected_transformation = Matrix4::<f32>::from_row_slice(&[
        angle_axis_rotation.angle.cos(), 0.0, angle_axis_rotation.angle.sin(), 0.0,
        0.0, 1.0, 0.0, 0.0,
        -angle_axis_rotation.angle.sin(), 0.0, angle_axis_rotation.angle.cos(), 0.0,
        0.0, 0.0, 0.0, 1.0                            ,
    ]);
    let errors = &( transformation - expected_transformation );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_exp_only_rotation_z() {
    let angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(0.0, 0.0, 1.0)).into_inner(),
        angle: 0.5
    };
    let twist = Twist::from_axis_angle_and_velocities(
        &angle_axis_rotation,
        &Vector3::new(0.0, 0.0, 0.0)
    );
    let algebra = twist.to_algebra();
    let transformation = algebra.exponential();
    let expected_transformation = Matrix4::<f32>::from_row_slice(&[
        angle_axis_rotation.angle.cos(), -angle_axis_rotation.angle.sin(), 0.0, 0.0,
        angle_axis_rotation.angle.sin(),  angle_axis_rotation.angle.cos(), 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0                            ,
    ]);
    let errors = &( transformation - expected_transformation );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_exp_rotation_and_axis() {
    let angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(0.0, 0.0, 1.0)).into_inner(),
        angle: PI
    };
    let twist = Twist::from_axis_angle_and_position_rotation(
        &angle_axis_rotation,
        &(Vector3::new(1.0, 0.0, 0.0))
    );
    let algebra = twist.to_algebra();
    let transformation = algebra.exponential();
    let expected_transformation = Matrix4::<f32>::from_row_slice(&[
        angle_axis_rotation.angle.cos(), -angle_axis_rotation.angle.sin(), 0.0,  2.0,
        angle_axis_rotation.angle.sin(),  angle_axis_rotation.angle.cos(), 0.0,  0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
    ]);
    let errors = &( transformation - expected_transformation );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_exp_rotation_z_and_translation() {
    let angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(0.0, 0.0, 1.0)).into_inner(),
        angle: PI
    };
    let twist = Twist::from_axis_angle_and_velocities(
        &angle_axis_rotation,
        &(Vector3::new(0.0, -1.0, 0.0) * PI)
    );
    let algebra = twist.to_algebra();
    let transformation = algebra.exponential();
    let expected_transformation = Matrix4::<f32>::from_row_slice(&[
        angle_axis_rotation.angle.cos(), -angle_axis_rotation.angle.sin(), 0.0,  2.0,
        angle_axis_rotation.angle.sin(),  angle_axis_rotation.angle.cos(), 0.0,  0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
    ]);
    let errors = &( transformation - expected_transformation );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_exp() {
    let algebra = Matrix4::<f32>::from_row_slice(&[
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, -1.57079632, 2.35619449,
        0.0, 1.57079632, 0.0, 2.35619449,
        0.0, 0.0, 0.0, 0.0,
    ]);
    let expected_transformation = Matrix4::<f32>::from_row_slice(&[
        1.0, 0.0,  0.0, 0.0,
        0.0, 0.0, -1.0, 0.0,
        0.0, 1.0,  0.0, 3.0,
        0.0, 0.0, 0.0,  1.0,
    ]);
    let transformation = algebra.exponential();
    let errors = &( transformation - expected_transformation );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_log() {
    let expected_algebra = Matrix4::<f32>::from_row_slice(&[
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, -1.57079632, 2.35619449,
        0.0, 1.57079632, 0.0, 2.35619449,
        0.0, 0.0, 0.0, 0.0,
    ]);
    let transformation = expected_algebra.exponential();
    let reconstructed_algebra = transformation.logarithm();
    let errors = &( reconstructed_algebra - expected_algebra );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_log_angle_axis_velocity() {
    let angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(1.0, 1.0, 1.0)).into_inner(),
        angle: PI / 3.0
    };
    let twist = Twist::from_axis_angle_and_velocities(
        &angle_axis_rotation,
        &(Vector3::new(3.0, -1.0, -4.0) * PI)
    );
    let expected_algebra = twist.to_algebra();
    let transformation = expected_algebra.exponential();
    let reconstructed_algebra = transformation.logarithm();
    let errors = &( reconstructed_algebra - expected_algebra );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_log_angle_axis_velocity_second() {
    let angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(2.0, 0.0, 1.0)).into_inner(),
        angle: PI * 1.0
    };
    let twist = Twist::from_axis_angle_and_velocities(
        &angle_axis_rotation,
        &(Vector3::new(0.0, -1.0, -4.0) * PI)
    );
    let expected_algebra = twist.to_algebra();
    let transformation = expected_algebra.exponential();
    let reconstructed_algebra = transformation.logarithm();
    let errors = &( reconstructed_algebra - expected_algebra );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_log_angle_axis_velocity_third() {
    let angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(1.0, 0.0, 0.0)).into_inner(),
        angle: PI * 2.0
    };
    let twist = Twist::from_axis_angle_and_velocities(
        &angle_axis_rotation,
        &(Vector3::new(0.0, 0.0, 0.0))
    );
    let expected_algebra = twist.to_algebra();
    let transformation = expected_algebra.exponential();
    let reconstructed_algebra = transformation.logarithm();
    let transformation_reconstructed = reconstructed_algebra.exponential();
    let transform = transformation_reconstructed.invert() * transformation;
    let errors = &( transform - Matrix4::<f32>::identity() );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_invert() {
    let angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(1.0, 3.0, 1.0)).into_inner(),
        angle: PI * 1.6
    };
    let twist = Twist::from_axis_angle_and_velocities(
        &angle_axis_rotation,
        &(Vector3::new(40.0, 10.0, -80.0))
    );
    let expected_algebra = twist.to_algebra();
    let transformation = expected_algebra.exponential();
    let transform = transformation.invert() * transformation;
    let errors = &( transform - Matrix4::<f32>::identity() );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_adjoint_is_equal_to_algebra_left_multiplied_and_inverse_right_multiplied() {
    let angle_axis_rotation = AxisAngleRotation {
        axis: Unit::new_normalize(Vector3::new(2.0, 1.0, 1.0)).into_inner(),
        angle: 1.6
    };
    let twist = Twist::from_axis_angle_and_velocities(
        &angle_axis_rotation,
        &(Vector3::new(10.0, -20.0, 80.0))
    );
    let algebra = twist.to_algebra();
    let transformation = algebra.exponential();
    let adjoint_transformation = transformation.to_adjoint();
    let twist_squared_computed = adjoint_transformation * twist;
    let expected_twist = (transformation * algebra * transformation.invert()).to_twist();
    let errors = &( twist_squared_computed - expected_twist );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_screw_chain() {
    let coordinates = &[0.0, -PI / 2., 0.0, 0.0, PI / 2.0, 0.0];
    let universal_robots_ur5 = UniversalRobotsUR5::from_default().unwrap();
    let transformation = universal_robots_ur5.to_transform(coordinates).unwrap();
    let expectd_transformation = Matrix4::<f32>::from_row_slice(&[
        0.0, -1.0, 0.0, 0.095,
        1.0,  0.0, 0.0, 0.109,
        0.0,  0.0, 1.0, 0.988,
        0.0,  0.0, 0.0, 1.0,
    ]);
    let errors = &( transformation - expectd_transformation );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}


#[test]
fn test_screw_chain_revolute_3_prismatic_1() {
    let coordinates = &[PI / 2.0, -PI, PI / 2.0, 10.0];
    let revolute_3_prismatic_1 = Revolute3Prismatic::from_default().unwrap();
    let transformation = revolute_3_prismatic_1.to_transform(coordinates).unwrap();
    let expectd_transformation = Matrix4::<f32>::from_row_slice(&[
            1.0, 0.0, 0.0,  0.0,
            0.0, 1.0, 0.0,  0.0,
            0.0, 0.0, 1.0, 10.0,
            0.0, 0.0, 0.0,  1.0,
        ]);
    let errors = &( transformation - expectd_transformation );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}

#[test]
fn test_jacobian_revolute_3_prismatic_1() {
    let coordinates = &[0.0, 0.0, 0.0, 0.0];
    let revolute_3_prismatic_1 = Revolute3Prismatic::from_default().unwrap();
    let jacobian = revolute_3_prismatic_1.to_space_jacobian(coordinates).unwrap();
    let expectd_jacobian = DMatrix::from_row_slice(6, coordinates.len(), &[
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            1.0, 1.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0,-1.0,-2.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
        ]);
    let errors = &( jacobian - expectd_jacobian );
    let error = errors.fold(0.0, |sum, element| sum + ( element * element ) );
    assert!(error < EPSILLON);
}