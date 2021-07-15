
pub mod algebraic_robots {

    use nalgebra::{Vector3, Vector6, Matrix4, U1, U3, U6, Matrix, SliceStorage, Matrix3, Matrix6, DMatrix, Unit};
    use std::cmp::{PartialOrd};
    use std::f64::consts::PI;

    pub static EPSILLON: f64 = 1e-14;
    pub type Vector3Slice<'a, T, RStride = U1, CStride = U6> =
        Matrix<T, U3, U1, SliceStorage<'a, T, U3, U1, RStride, CStride>>;

    #[derive(Debug, PartialEq, PartialOrd)]
    pub struct AxisAngleRotation {
        pub axis: Vector3<f64>,
        pub angle: f64
    }

    pub type Twist  = Vector6<f64>;
    pub type Screw  = Vector6<f64>;
    pub type Wrench = Vector6<f64>;
    pub type ProjectiveAlgebraRep = Matrix4<f64>;
    pub type ProjectiveGroupRep   = Matrix4<f64>;
    pub type ProjectiveAdjointRep = Matrix6<f64>;

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
        fn angular(&self) -> Vector3Slice<f64>;
        fn linear(&self)  -> Vector3Slice<f64>;
        fn from_angular_linear(angular : Vector3<f64>, linear : Vector3<f64>) -> Vector6<f64>;
    }

    pub trait TwistLike {
        fn from_axis_angle_and_position_rotation(
            axis_angle_rotation: &AxisAngleRotation,
            axis_point: &Vector3<f64>) -> Twist;
        fn from_axis_angle_and_velocities(
                axis_angle_rotation: &AxisAngleRotation,
                linear_velocity: &Vector3<f64>) -> Twist;
        fn to_axis_angle_rotation(&self) -> AxisAngleRotation;
        fn to_algebra(&self) -> ProjectiveAlgebraRep;
        fn to_screw(&self) -> Screw;
    }

    pub trait WrenchLike {
        fn from_point_and_force(point: Vector3<f64>, force: Vector3<f64>) -> Wrench;
        fn from_moment_and_force(moment: Vector3<f64>, force: Vector3<f64>) -> Wrench;
    }


    impl ScrewChain {

        pub fn to_transform(&self, coordinates : &[f64]) -> Option<ProjectiveGroupRep> {
            if coordinates.len() == self.screws.len() {
                let transform_screws = self.screws.iter().zip(coordinates.iter()).
                    fold( Matrix4::identity(),
                        | current_transform : Matrix4<f64>, (&screw, &coordinate) |
                        current_transform * ( (screw * coordinate).to_algebra().exponential() ));
                return Option::Some(transform_screws * self.end_effector_at_initial_position )
            } else {
                return Option::None
            }
        }

        fn add_transformation(transforms : Vec<Matrix4<f64>>,
                screw_coordinate: (&Vector6<f64>, &f64) ) -> Vec<Matrix4<f64>> {
            let (screw, coordinate) = screw_coordinate;
            let screw_transformation : Matrix4<f64> = (*screw * *coordinate).to_algebra().exponential();
            let new_transformation : Matrix4<f64> = transforms.last().unwrap() * screw_transformation;
            let new_transformations = [transforms, vec![new_transformation]].concat();
            new_transformations
        }

        fn add_jacobinan_column(jacobian_columns : Vec<Vector6<f64>>,
                screw_and_incremental_transformation: (&Vector6<f64>, &Matrix4<f64>) ) -> Vec<Vector6<f64>> {
            let (screw,_transformation) = screw_and_incremental_transformation;
            let adjoint_transformation : Matrix6<f64> = _transformation.to_adjoint();
            let jacobian_column : Vector6<f64> = adjoint_transformation *  screw;
            let new_jacobian_columns = [jacobian_columns, vec![jacobian_column]].concat();
            new_jacobian_columns
        }

        pub fn to_space_jacobian(&self, coordinates : &[f64]) -> Option<DMatrix<f64> > {
            if coordinates.len() == self.screws.len() {
                let screws_and_coordinates = self.screws.iter().zip(coordinates.iter());
                let identity : Vec<Matrix4<f64>> = vec![Matrix4::identity()];
                let incremental_transform_screws : Vec<Matrix4<f64>> = screws_and_coordinates.fold(
                    identity, ScrewChain::add_transformation);
                let screws_and_incremental_transformations = self.screws.iter().zip(
                    incremental_transform_screws.iter());
                let jacobian_columns : &Vec<Vector6<f64>> = &screws_and_incremental_transformations.fold(
                    vec![], ScrewChain::add_jacobinan_column);
                let raw_iter = jacobian_columns.into_iter().flat_map(|vector| vector.into_iter()).
                    map(|x| *x).collect::<Vec<_>>();
                let jacobian_matrix : DMatrix<f64> = DMatrix::from_vec(6, coordinates.len(), raw_iter);
                return Option::Some(jacobian_matrix)
            } else {
                return Option::None
            }
        }

    }

    impl WrenchLike for Wrench {

        fn from_point_and_force(point: Vector3<f64>, force: Vector3<f64>) -> Wrench {
            Vector6::<f64>::from_row_slice( &[ point.cross(&force).as_slice(), force.as_slice() ].concat() )
        }

        fn from_moment_and_force(moment: Vector3<f64>, force: Vector3<f64>) -> Wrench {
            Vector6::<f64>::from_row_slice( &[ moment.as_slice(), force.as_slice() ].concat() )
        }

    }

    impl GeneralizedCoordinates for Vector6<f64> {

        fn from_angular_linear(angular : Vector3<f64>, linear : Vector3<f64>) -> Vector6<f64> {
            Vector6::<f64>::from_row_slice(&[angular.as_slice(), linear.as_slice()].concat())
        }

        fn angular(&self) -> Vector3Slice<f64> {
            self.fixed_slice::<3, 1>(0, 0)
        }

        fn linear(&self) -> Vector3Slice<f64> {
            self.fixed_slice::<3, 1>(3, 0)
        }

    }

    impl SE3Group for ProjectiveGroupRep {

        fn to_adjoint(&self) -> ProjectiveAdjointRep {
            let so3 = self.fixed_slice::<3, 3>(0, 0);
            let translation = self.fixed_slice::<3, 1>(0, 3);
            let skew_translation =  Matrix3::<f64>::from_row_slice(&[
                           0.0, -translation[2],  translation[1],
                 translation[2],            0.0, -translation[0],
                -translation[1],  translation[0],            0.0
                ]);
            let t3 = skew_translation * so3;
            Matrix6::<f64>::from_row_slice(&[
                *so3.index((0, 0)), *so3.index((0, 1)), *so3.index((0, 2)),                0.0,                0.0,                0.0,
                *so3.index((1, 0)), *so3.index((1, 1)), *so3.index((1, 2)),                0.0,                0.0,                0.0,
                *so3.index((2, 0)), *so3.index((2, 1)), *so3.index((2, 2)),                0.0,                0.0,                0.0,
                 *t3.index((0, 0)),  *t3.index((0, 1)),  *t3.index((0, 2)), *so3.index((0, 0)), *so3.index((0, 1)), *so3.index((0, 2)),
                 *t3.index((1, 0)),  *t3.index((1, 1)),  *t3.index((1, 2)), *so3.index((1, 0)), *so3.index((1, 1)), *so3.index((1, 2)),
                 *t3.index((2, 0)),  *t3.index((2, 1)),  *t3.index((2, 2)), *so3.index((2, 0)), *so3.index((2, 1)), *so3.index((2, 2)),
            ])
        }

        fn logarithm(&self) -> ProjectiveAlgebraRep {
            let so3 = self.fixed_slice::<3, 3>(0, 0);
            let acosinput =  {
                let acosinput_raw = (so3.trace() - 1.0) / 2.0;
                if (acosinput_raw + 1.0).abs() < EPSILLON {
                    -1.0
                } else if (acosinput_raw - 1.0).abs() < EPSILLON {
                    1.0
                } else {
                    acosinput_raw
                }
            };
            let theta = acosinput.acos();
            let so3_alg = {
                if acosinput >= 1.0 {
                    Matrix3::<f64>::zeros()
                } else if acosinput <= -1.0 {
                    let omg = {
                        if  ( 1.0 + so3.index((2, 2)) ).abs() > EPSILLON {
                            ( 1.0 / ( 2.0 * ( 1.0 + so3.index((2, 2)))).sqrt() )
                                * Vector3::new( *so3.index((0, 2)), *so3.index((1, 2)), 1.0 + (*so3.index((2, 2))) )
                        } else if (1.0 + so3.index((1, 1))).abs() > EPSILLON {
                            ( 1.0 / ( 2.0 * (1.0 + *so3.index((1, 1)) )).sqrt() )
                                * Vector3::new( *so3.index((0, 1)), 1.0  + *so3.index((1, 1)), *so3.index((2,1)) )
                        } else {
                            (1.0 / ( 2.0 * ( 1.0 + *so3.index((0, 0))) ).sqrt() )
                                * Vector3::new( 1.0 + *so3.index((0, 0)), *so3.index((1, 0)), *so3.index((2, 0)))
                        }
                    };
                    Matrix3::<f64>::from_row_slice(&[
                            0.0, -omg[2],  omg[1],
                         omg[2],     0.0, -omg[0],
                        -omg[1],  omg[0],     0.0
                    ]) * PI
                } else {
                    ( ( theta / 2.0 ) / theta.sin() ) * ( so3 - so3.transpose() )
                }
            };
            if so3_alg == Matrix3::<f64>::zeros() {
                Matrix4::<f64>::from_row_slice(&[
                    0.0, 0.0, 0.0, *self.index((0, 3)),
                    0.0, 0.0, 0.0, *self.index((1, 3)),
                    0.0, 0.0, 0.0, *self.index((2, 3)),
                    0.0, 0.0, 0.0,                 0.0,
                ])
            } else {
                let so3_algebra_square = so3_alg * so3_alg;
                let t3_alg = (
                                    Matrix3::<f64>::identity() -
                                    ( so3_alg / 2.0 ) +
                                    (
                                        ( 1.0 / theta ) -
                                        ( ( 1.0 / ( theta / 2.0 ).tan() ) / 2.0 )
                                    ) * ( so3_algebra_square / theta )
                                ) * self.fixed_slice::<3,1>(0, 3);
                Matrix4::<f64>::from_row_slice(&[
                    *so3_alg.index((0, 0)), *so3_alg.index((0, 1)), *so3_alg.index((0, 2)), *t3_alg.index((0, 0)),
                    *so3_alg.index((1, 0)), *so3_alg.index((1, 1)), *so3_alg.index((1, 2)), *t3_alg.index((1, 0)),
                    *so3_alg.index((2, 0)), *so3_alg.index((2, 1)), *so3_alg.index((2, 2)), *t3_alg.index((2, 0)),
                                       0.0,                    0.0,                    0.0,                   0.0,
                ])
            }
        }

        fn invert(&self) -> ProjectiveGroupRep {
            let rotation = self.fixed_slice::<3, 3>(0, 0);
            let invert_so3 = rotation.transpose();
            let translation = self.fixed_slice::<3, 1>(0, 3);
            let inv_t3 = - invert_so3 * translation;
            Matrix4::<f64>::from_row_slice(&[
                *invert_so3.index((0, 0)), *invert_so3.index((0, 1)), *invert_so3.index((0, 2)), *inv_t3.index((0, 0)),
                *invert_so3.index((1, 0)), *invert_so3.index((1, 1)), *invert_so3.index((1, 2)), *inv_t3.index((1, 0)),
                *invert_so3.index((2, 0)), *invert_so3.index((2, 1)), *invert_so3.index((2, 2)), *inv_t3.index((2, 0)),
                                      0.0,                       0.0,                       0.0,                   1.0,
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
                let so3 = Matrix3::<f64>::identity() + ( angle_sin * omgmat) +
                     ( ( 1.0 - angle_cos ) * omgmat_square );
                let t3_group = (
                                ( Matrix3::<f64>::identity() * axis_angle_rotation.angle ) +
                                ( ( 1.0 - angle_cos ) * omgmat ) +
                                ( ( axis_angle_rotation.angle - angle_sin ) * omgmat_square )
                             ) * ( self.fixed_slice::<3, 1>(0, 3) / axis_angle_rotation.angle );
                Matrix4::<f64>::from_row_slice(&[
                    *so3.index((0, 0)), *so3.index((0, 1)), *so3.index((0, 2)), *t3_group.index((0, 0)),
                    *so3.index((1, 0)), *so3.index((1, 1)), *so3.index((1, 2)), *t3_group.index((1, 0)),
                    *so3.index((2, 0)), *so3.index((2, 1)), *so3.index((2, 2)), *t3_group.index((2, 0)),
                                   0.0,                0.0,                0.0,                     1.0,
                ])
            } else {
                let linear_velocity = twist.linear();
                Matrix4::<f64>::from_row_slice(&[
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
            axis_point: &Vector3<f64>) -> Twist {
            let angular_velocity = axis_angle_rotation.angle * axis_angle_rotation.axis;
            let linear_velocity = - angular_velocity.cross(axis_point);
            Twist::from_angular_linear(
                angular_velocity,
                linear_velocity
            )
        }

        fn from_axis_angle_and_velocities(
            axis_angle_rotation: &AxisAngleRotation,
            linear_velocity: &Vector3<f64>) -> Twist {
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
                    Vector3::<f64>::from_row_slice(angular_velocity.as_slice())).into_inner(),
                angle: angle
            }
        }

        fn to_algebra(&self) -> ProjectiveAlgebraRep {
            let angular_velocity = self.angular();
            let linear_velocity = self.linear();
            Matrix4::<f64>::from_row_slice(&[
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

}

pub mod screw_chains;

#[cfg(test)]
mod test;
