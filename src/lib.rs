
pub mod algebraic_robots {

    use nalgebra::{Vector3, Vector6, Matrix4, U1, U3, U6, Matrix, SliceStorage, Matrix3, Matrix6, Unit};
    use std::cmp::{PartialOrd};
    use std::f32::consts::PI;

    pub static EPSILLON: f32 = 1e-6;
    pub type Vector3Slice<'a, T, RStride = U1, CStride = U6> = Matrix<T, U3, U1, SliceStorage<'a, T, U3, U1, RStride, CStride>>;

    #[derive(Debug, PartialEq, PartialOrd)]
    pub struct AxisAngleRotation {
        pub axis: Vector3<f32>,
        pub angle: f32
    }

    pub type Twist  = Vector6<f32>;
    pub type Screw  = Vector6<f32>;
    pub type Wrench = Vector6<f32>;
    pub type ProjectiveAlgebraRep = Matrix4<f32>;
    pub type ProjectiveGroupRep   = Matrix4<f32>;
    pub type ProjectiveAdjointRep = Matrix6<f32>;

    pub struct ScrewChain {
        pub screws: Vec<Screw>,
        pub end_effector_at_initial_position: ProjectiveGroupRep,
    }

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

    pub trait WrenchLike {
        fn from_point_and_force(point: Vector3<f32>, force: Vector3<f32>) -> Wrench;
        fn from_moment_and_force(moment: Vector3<f32>, force: Vector3<f32>) -> Wrench;
    }


    impl ScrewChain {

        pub fn to_transform(&self, coordinates : &[f32]) -> Option<ProjectiveGroupRep> {
            if coordinates.len() == self.screws.len() {
                let transform_screws = self.screws.iter().zip(coordinates.iter()).
                    fold( Matrix4::identity(),
                        | current_transform : Matrix4<f32>, (&screw, &coordinate) |
                        current_transform * ( (screw * coordinate).to_algebra().exponential() ));
                return Option::Some(transform_screws * self.end_effector_at_initial_position )
            } else {
                return Option::None
            }
        }

        pub fn to_space_jacobian(&self, coordinates : &[f32]) -> Option<ProjectiveAdjointRep> {
            if coordinates.len() == self.screws.len() {
                let transform_screws = self.screws.iter().rev().zip(coordinates.iter().rev()).
                    fold( Matrix6::identity(),
                        | current_transform : Matrix6<f32>, (&screw, &coordinate) |
                        current_transform * ( (screw * coordinate).to_algebra().exponential().to_adjoint() ));
                return Option::Some(transform_screws)
            } else {
                return Option::None
            }

        }

        pub fn to_body_jacobian_TODO(&self, coordinates : &[f32]) -> Option<ProjectiveAdjointRep> {
            if coordinates.len() == self.screws.len() {
                let transform_screws = self.screws.iter().rev().zip(coordinates.iter().rev()).
                    fold( Matrix6::identity(),
                        | current_transform : Matrix6<f32>, (&screw, &coordinate) |
                        current_transform * ( (screw * coordinate).to_algebra().exponential().to_adjoint() ));
                return Option::Some(transform_screws)
            } else {
                return Option::None
            }

        }

    }

    impl WrenchLike for Wrench {

        fn from_point_and_force(point: Vector3<f32>, force: Vector3<f32>) -> Wrench {
            Vector6::<f32>::from_row_slice( &[ point.cross(&force).as_slice(), force.as_slice() ].concat() )
        }

        fn from_moment_and_force(moment: Vector3<f32>, force: Vector3<f32>) -> Wrench {
            Vector6::<f32>::from_row_slice( &[ moment.as_slice(), force.as_slice() ].concat() )
        }

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
                    0.0, 0.0, 0.0,                 0.0,
                ])
            } else {
                let so3_algebra_square = so3_algebra * so3_algebra;
                let t3_algebra = (
                                    Matrix3::<f32>::identity() -
                                    ( so3_algebra / 2.0 ) +
                                    (
                                        ( 1.0 / theta ) -
                                        ( ( 1.0 / ( theta / 2.0 ).tan() ) / 2.0 )
                                    ) * ( so3_algebra_square / theta )
                                ) * self.fixed_slice::<3,1>(0, 3);
                Matrix4::<f32>::from_row_slice(&[
                    *so3_algebra.index((0, 0)), *so3_algebra.index((0, 1)), *so3_algebra.index((0, 2)), *t3_algebra.index((0, 0)),
                    *so3_algebra.index((1, 0)), *so3_algebra.index((1, 1)), *so3_algebra.index((1, 2)), *t3_algebra.index((1, 0)),
                    *so3_algebra.index((2, 0)), *so3_algebra.index((2, 1)), *so3_algebra.index((2, 2)), *t3_algebra.index((2, 0)),
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
                                         0.0,                      0.0,                      0.0,                     1.0,
                ])
            } else {
                let linear_velocity = twist.linear();
                Matrix4::<f32>::from_row_slice(&[
                        1.0, 0.0, 0.0, linear_velocity[0],
                        0.0, 1.0, 0.0, linear_velocity[1],
                        0.0, 0.0, 1.0, linear_velocity[2],
                        0.0, 0.0, 0.0,                1.0,
                    ])

            }

        }

        fn to_twist(&self) -> Twist {
            Twist::from_angular_linear(
                    Vector3::new( *self.index((2, 1)), *self.index((0, 2)), *self.index((1, 0)) ),
                    Vector3::new( *self.index((0, 3)), *self.index((1, 3)), *self.index((2, 3)) )
                )
        }
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
                                     0.0, -angular_velocity[2],  angular_velocity[1], linear_velocity[0],
                     angular_velocity[2],                  0.0, -angular_velocity[0], linear_velocity[1],
                    -angular_velocity[1],  angular_velocity[0],                  0.0, linear_velocity[2],
                                     0.0,                  0.0,                  0.0,                0.0,
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
                let linear_velocity_norm = linear_velocity.norm();
                if linear_velocity_norm > EPSILLON {
                    Screw::from_angular_linear(
                        Vector3::new(0.0, 0.0, 0.0),
                        linear_velocity / linear_velocity_norm
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

    pub mod screw_chains {

        use super::*;

        pub mod universal_robot_ur5 {

            use super::*;

            pub fn create() -> ScrewChain {
                let h1 =  89.0 / 1000.0;
                let h2 =  95.0 / 1000.0;
                let l1 = 425.0 / 1000.0;
                let l2 = 392.0 / 1000.0;
                let w1 = 109.0 / 1000.0;
                let w2 =  82.0 / 1000.0;
                from_parameters(h1, h2, l1, l2, w1, w2)
            }

            pub fn from_parameters(h1: f32, h2: f32, l1: f32, l2: f32, w1: f32, w2: f32) -> ScrewChain {
                ScrewChain {
                    screws: vec![
                            Screw::from_angular_linear(  Vector3::new(0.0, 0.0, 1.0), Vector3::new(0.0, 0.0, 0.0)),
                            Screw::from_angular_linear(  Vector3::new(0.0, 1.0, 0.0), Vector3::new(-h1, 0.0, 0.0)),
                            Screw::from_angular_linear(  Vector3::new(0.0, 1.0, 0.0), Vector3::new(-h1, 0.0, l1)),
                            Screw::from_angular_linear(  Vector3::new(0.0, 1.0, 0.0), Vector3::new(-h1, 0.0, l1 + l2)),
                            Screw::from_angular_linear( Vector3::new(0.0, 0.0, -1.0), Vector3::new(-w1, l1 + l2, 0.0)),
                            Screw::from_angular_linear(  Vector3::new(0.0, 1.0, 0.0), Vector3::new(h2 - h1, 0.0, l1 + l2))
                        ],
                    end_effector_at_initial_position: Matrix4::<f32>::from_row_slice(&[
                        -1.0, 0.0, 0.0, l1 + l2,
                         0.0, 0.0, 1.0, w1 + w2,
                         0.0, 1.0, 0.0, h1 - h2,
                         0.0, 0.0, 0.0,     1.0,
                    ])
                }
            }

        }

    }
}

#[cfg(test)]
mod test;
