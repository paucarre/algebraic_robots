
pub mod algebraic_robots {

    use nalgebra::{Vector3, Vector6, Matrix4, U1, U3, U6, VectorSlice3, VectorSlice6, Matrix, SliceStorage, Matrix3, Matrix6, Unit};
    use std::cmp::{PartialOrd};
    use std::f32::consts::PI;

    pub static EPSILLON: f32 = 1e-6;
    pub type Vector3Slice<'a, T, RStride = U1, CStride = U6> = Matrix<T, U3, U1, SliceStorage<'a, T, U3, U1, RStride, CStride>>;

    #[derive(Debug, PartialEq, PartialOrd)]
    pub struct AxisAngleRotation {
        pub axis: Vector3<f32>,
        pub angle: f32
    }

    pub type Twist = Vector6<f32>;
    pub type Screw = Vector6<f32>;
    pub type ProjectiveAlgebraRep = Matrix4<f32>;
    pub type ProjectiveGroupRep = Matrix4<f32>;
    pub type ProjectiveAdjointRep = Matrix6<f32>;

    pub trait SE3Algebra {
        fn to_twist(&self) -> Twist;
        fn exponential(&self) -> ProjectiveAlgebraRep;
    }

    pub trait SE3Group {
        fn invert(&self) -> ProjectiveGroupRep;
        fn logarithm(&self) -> ProjectiveAlgebraRep;
        fn to_adjoint(&self) -> ProjectiveAdjointRep;
    }

    pub trait GeneralizedCoordinates {
        fn angular(&self) -> Vector3Slice<f32>;
        fn linear(&self)  -> Vector3Slice<f32>;
        fn from_angular_linear(angular : Vector3<f32>, linear : Vector3<f32>) -> Vector6<f32>;
    }

    impl GeneralizedCoordinates for Vector6<f32> {

        fn from_angular_linear(angular : Vector3<f32>, linear : Vector3<f32>) -> Vector6<f32> {
            Vector6::<f32>::from_row_slice(&[angular.as_slice(), linear.as_slice()].concat())
        }

        fn angular(&self) -> Vector3Slice<f32> {
            self.fixed_slice::<3, 1>(0, 0)
        }

        fn linear(&self) -> Vector3Slice<f32> {
            self.fixed_slice::<3, 1>(3, 0)
        }

    }

    impl SE3Group for ProjectiveGroupRep {

        fn to_adjoint(&self) -> ProjectiveAdjointRep {
            let rotation = self.fixed_slice::<3, 3>(0, 0);
            let translation = self.fixed_slice::<3, 1>(0, 3);
            let skew_translation =  Matrix3::<f32>::from_row_slice(&[
                           0.0, -translation[2],  translation[1],
                 translation[2],            0.0, -translation[0],
                -translation[1],  translation[0],            0.0
                ]);
            let bottom_left = skew_translation * rotation;
            Matrix6::<f32>::from_row_slice(&[
                   *rotation.index((0, 0)),    *rotation.index((0, 1)),    *rotation.index((0, 2)),                    0.0,                     0.0,                      0.0,
                   *rotation.index((1, 0)),    *rotation.index((1, 1)),    *rotation.index((1, 2)),                    0.0,                     0.0,                      0.0,
                   *rotation.index((2, 0)),    *rotation.index((2, 1)),    *rotation.index((2, 2)),                    0.0,                     0.0,                      0.0,
                *bottom_left.index((0, 0)), *bottom_left.index((0, 1)), *bottom_left.index((0, 2)), *rotation.index((0, 0)), *rotation.index((0, 1)), *rotation.index((0, 2)),
                *bottom_left.index((1, 0)), *bottom_left.index((1, 1)), *bottom_left.index((1, 2)), *rotation.index((1, 0)), *rotation.index((1, 1)), *rotation.index((1, 2)),
                *bottom_left.index((2, 0)), *bottom_left.index((2, 1)), *bottom_left.index((2, 2)), *rotation.index((2, 0)), *rotation.index((2, 1)), *rotation.index((2, 2)),
            ])
        }

        fn logarithm(&self) -> ProjectiveAlgebraRep {
            let rotation = self.fixed_slice::<3, 3>(0, 0);
            let acosinput =  {
                let acosinput_raw = (rotation.trace() - 1.0) / 2.0;
                if (acosinput_raw + 1.0).abs() < EPSILLON {
                    -1.0
                } else if (acosinput_raw - 1.0).abs() < EPSILLON {
                    1.0
                } else {
                    acosinput_raw
                }
            };
            let theta = acosinput.acos();
            let so3_algebra = {
                if acosinput >= 1.0 {
                    Matrix3::<f32>::zeros()
                } else if acosinput <= -1.0 {
                    let omg = {
                        if  ( 1.0 + rotation.index((2, 2)) ).abs() > EPSILLON {
                            ( 1.0 / ( 2.0 * ( 1.0 + rotation.index((2, 2)))).sqrt() )
                                * Vector3::new( *rotation.index((0, 2)), *rotation.index((1, 2)), 1.0 + (*rotation.index((2, 2))) )
                        } else if (1.0 + rotation.index((1, 1))).abs() > EPSILLON {
                            ( 1.0 / ( 2.0 * (1.0 + *rotation.index((1, 1)) )).sqrt() )
                                * Vector3::new( *rotation.index((0, 1)), 1.0  + *rotation.index((1, 1)), *rotation.index((2,1)) )
                        } else {
                            (1.0 / ( 2.0 * ( 1.0 + *rotation.index((0, 0))) ).sqrt() )
                                * Vector3::new( 1.0 + *rotation.index((0, 0)), *rotation.index((1, 0)), *rotation.index((2, 0)))
                        }
                    };
                    Matrix3::<f32>::from_row_slice(&[
                            0.0, -omg[2],  omg[1],
                        omg[2],     0.0, -omg[0],
                        -omg[1],  omg[0],     0.0
                    ]) * PI
                } else {
                    ( ( theta / 2.0 ) / theta.sin() ) * ( rotation - rotation.transpose() )
                }
            };
            if so3_algebra == Matrix3::<f32>::zeros() {
                Matrix4::<f32>::from_row_slice(&[
                    0.0, 0.0, 0.0, *self.index((0, 3)),
                    0.0, 0.0, 0.0, *self.index((1, 3)),
                    0.0, 0.0, 0.0, *self.index((2, 3)),
                    0.0, 0.0, 0.0,                  0.0,
                ])
            } else {
                let so3_algebra_square = so3_algebra * so3_algebra;
                let velocities = (
                                    Matrix3::<f32>::identity() -
                                    ( so3_algebra / 2.0 ) +
                                    (
                                        ( 1.0 / theta ) -
                                        ( ( 1.0 / ( theta / 2.0 ).tan() ) / 2.0 )
                                    ) * ( so3_algebra_square / theta )
                                ) * self.fixed_slice::<3,1>(0, 3);
                Matrix4::<f32>::from_row_slice(&[
                    *so3_algebra.index((0, 0)), *so3_algebra.index((0, 1)), *so3_algebra.index((0, 2)), *velocities.index((0, 0)),
                    *so3_algebra.index((1, 0)), *so3_algebra.index((1, 1)), *so3_algebra.index((1, 2)), *velocities.index((1, 0)),
                    *so3_algebra.index((2, 0)), *so3_algebra.index((2, 1)), *so3_algebra.index((2, 2)), *velocities.index((2, 0)),
                                           0.0,                        0.0,                        0.0,                       0.0,
                ])
            }
        }

        fn invert(&self) -> ProjectiveGroupRep {
            let rotation = self.fixed_slice::<3, 3>(0, 0);
            let inverted_rotation = rotation.transpose();
            let translation = self.fixed_slice::<3, 1>(0, 3);
            let interted_translation = - inverted_rotation * translation;
            Matrix4::<f32>::from_row_slice(&[
                *inverted_rotation.index((0, 0)), *inverted_rotation.index((0, 1)), *inverted_rotation.index((0, 2)), *interted_translation.index((0, 0)),
                *inverted_rotation.index((1, 0)), *inverted_rotation.index((1, 1)), *inverted_rotation.index((1, 2)), *interted_translation.index((1, 0)),
                *inverted_rotation.index((2, 0)), *inverted_rotation.index((2, 1)), *inverted_rotation.index((2, 2)), *interted_translation.index((2, 0)),
                                             0.0,                              0.0,                              0.0,                                 1.0,
            ])
        }

    }



    impl SE3Algebra for ProjectiveAlgebraRep {

        fn exponential(&self) -> ProjectiveGroupRep {
            let twist = self.to_twist();
            let angular_velocity = twist.angular();
            if angular_velocity.norm() > EPSILLON {
                let axis_angle_rotation = twist.to_axis_angle_rotation();
                let rotation_algebra = self.fixed_slice::<3, 3>(0, 0);
                let omgmat = rotation_algebra / axis_angle_rotation.angle;
                let omgmat_square = omgmat * omgmat;
                let angle_sin = axis_angle_rotation.angle.sin();
                let angle_cos = axis_angle_rotation.angle.cos();
                let so3_group = Matrix3::<f32>::identity() + ( angle_sin * omgmat) +
                     ( ( 1.0 - angle_cos ) * omgmat_square );
                let t3_group = (
                                ( Matrix3::<f32>::identity() * axis_angle_rotation.angle ) +
                                ( ( 1.0 - angle_cos ) * omgmat ) +
                                ( ( axis_angle_rotation.angle - angle_sin ) * omgmat_square )
                             ) * ( self.fixed_slice::<3, 1>(0, 3) / axis_angle_rotation.angle );
                Matrix4::<f32>::from_row_slice(&[
                    *so3_group.index((0, 0)), *so3_group.index((0, 1)), *so3_group.index((0, 2)), *t3_group.index((0, 0)),
                    *so3_group.index((1, 0)), *so3_group.index((1, 1)), *so3_group.index((1, 2)), *t3_group.index((1, 0)),
                    *so3_group.index((2, 0)), *so3_group.index((2, 1)), *so3_group.index((2, 2)), *t3_group.index((2, 0)),
                    0.0, 0.0, 0.0, 1.0                            ,
                ])
            } else {
                let linear_velocity = twist.linear();
                Matrix4::<f32>::from_row_slice(&[
                        1.0, 0.0, 0.0, linear_velocity[0],
                        0.0, 1.0, 0.0, linear_velocity[1],
                        0.0, 0.0, 1.0, linear_velocity[2],
                        0.0, 0.0, 0.0, 1.0                            ,
                    ])

            }

        }

        fn to_twist(&self) -> Twist {
            Twist::from_angular_linear(
                Vector3::new(
                    *self.index((2, 1)),
                    *self.index((0, 2)),
                    *self.index((1, 0))),
                Vector3::new(
                    *self.index((0, 3)),
                    *self.index((1, 3)),
                    *self.index((2, 3)))
                )
        }
    }

    pub trait TwistLike {
        fn from_axis_angle_and_position_rotation(
            axis_angle_rotation: &AxisAngleRotation,
            axis_point: &Vector3<f32>) -> Twist;
        fn from_axis_angle_and_velocities(
                axis_angle_rotation: &AxisAngleRotation,
                linear_velocity: &Vector3<f32>) -> Twist;
        fn to_axis_angle_rotation(&self) -> AxisAngleRotation;
        fn to_algebra(&self) -> ProjectiveAlgebraRep;
        fn to_screw(&self) -> Screw;
    }

    impl TwistLike for Twist {

        fn from_axis_angle_and_position_rotation(
            axis_angle_rotation: &AxisAngleRotation,
            axis_point: &Vector3<f32>) -> Twist {
            let angular_velocity = axis_angle_rotation.angle * axis_angle_rotation.axis;
            let linear_velocity = - angular_velocity.cross(axis_point);
            Twist::from_angular_linear(
                angular_velocity,
                linear_velocity
            )
        }

        fn from_axis_angle_and_velocities(
            axis_angle_rotation: &AxisAngleRotation,
            linear_velocity: &Vector3<f32>) -> Twist {
            Twist::from_angular_linear(
                axis_angle_rotation.angle * axis_angle_rotation.axis,
                *linear_velocity
            )
        }

        fn to_axis_angle_rotation(&self) -> AxisAngleRotation {
            let angular_velocity = self.angular();
            let angle = angular_velocity.norm();
            AxisAngleRotation {
                axis: Unit::new_normalize(
                    Vector3::<f32>::from_row_slice(angular_velocity.as_slice())).into_inner(),
                angle: angle
            }
        }

        fn to_algebra(&self) -> ProjectiveAlgebraRep {
            let angular_velocity = self.angular();
            let linear_velocity = self.linear();
            Matrix4::<f32>::from_row_slice(&[
                    0.0                      , -angular_velocity[2],  angular_velocity[1], linear_velocity[0],
                    angular_velocity[2] , 0.0                      , -angular_velocity[0], linear_velocity[1],
                    -angular_velocity[1], angular_velocity[0] ,  0.0                     , linear_velocity[2],
                    0.0                      , 0.0                      ,  0.0                     , 0.0                           ,
                ])
        }

        fn to_screw(&self) -> Screw {
            let angular_velocity = self.angular();
            let linear_velocity = self.linear();
            let angle = angular_velocity.norm();
            if angle > EPSILLON {
                Screw::from_angular_linear(
                     angular_velocity / angle,
                     linear_velocity / angle
                )
            } else {
                let translation = linear_velocity.norm();
                if translation > EPSILLON {
                    Screw::from_angular_linear(
                        Vector3::new(0.0, 0.0, 0.0),
                        linear_velocity / translation
                    )
                } else {
                    Screw::from_angular_linear(
                        Vector3::new(0.0, 0.0, 0.0),
                        Vector3::new(0.0, 0.0, 0.0)
                    )
                }
            }
        }

    }

    #[cfg(test)]
    mod tests {

        use super::*;
        use { EPSILLON, AxisAngleRotation, Twist, SE3Algebra };
        use nalgebra::{ Unit, Vector3 };
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
                angle, 0.0, -angle, 30.0,
                -angle, angle,  0.0, 40.0,
                0.0, 0.0,  0.0, 0.0,
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
                1.0, 0.0, 0.0, 20.0,
                0.0, 1.0, 0.0, -30.0,
                0.0, 0.0, 1.0, 40.0,
                0.0, 0.0, 0.0, 1.0                            ,
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
                angle_axis_rotation.angle.sin(), angle_axis_rotation.angle.cos(), 0.0, 0.0,
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
        fn test_adjoint() {
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

    }
}


