use crate::algebraic_robots::{Screw, GeneralizedCoordinates, ScrewChain};
use nalgebra::{ Vector3, Matrix4};
use std::collections::HashMap;

pub trait ScrewChainLike {
    fn from_default() -> Option<ScrewChain>;
    fn from_parameters(topological_parameters: &HashMap<&str, f64>) -> Option<ScrewChain>;
    fn parameters() -> Vec<String>;
    fn param_map_to_param_vect(parameters_map: &HashMap<&str, f64>) -> Option<Vec<f64>> {
        let params = Self::parameters();
        let optional_values = params.iter().map(|parameter| parameters_map.get(&parameter.as_str()));
        let empty_vec = Some(vec![]);
        let params_opt = optional_values.fold(empty_vec, |opt_vect_params, parameter| match opt_vect_params {
                Some(vect_params) => match parameter {
                    Some(param) => Some([vect_params, vec![*param]].concat()),
                    None => None
                },
                None => None
        });
        params_opt
    }
}

pub struct UniversalRobotsUR5;
pub struct Revolute3PrismaticRevolute3;
pub struct Revolute3Prismatic;



impl ScrewChainLike for UniversalRobotsUR5 {

    fn from_default() -> Option<ScrewChain> {
        let default_parameters: HashMap<&str, f64> = [
            ("h1",  89.0 / 1000.0),
            ("h2",  95.0 / 1000.0),
            ("l1", 425.0 / 1000.0),
            ("l2", 392.0 / 1000.0),
            ("w1", 109.0 / 1000.0),
            ("w2",  82.0 / 1000.0)
            ].iter().cloned().collect();
        Self::from_parameters(&default_parameters)
    }

    fn parameters() -> Vec<String> {
        return vec![String::from("h1"), String::from("h2"),
                    String::from("l1"), String::from("l2"),
                    String::from("w1"), String::from("w2")];
    }

    fn from_parameters(parameters_map: &HashMap<&str, f64>) -> Option<ScrewChain> {
        let params_opt = Self::param_map_to_param_vect(parameters_map);
        params_opt.map(|params|
            ScrewChain {
                screws: vec![
                        Screw::from_angular_linear(  Vector3::new(0.0, 0.0,  1.0), Vector3::new(0.0, 0.0, 0.0)),
                        Screw::from_angular_linear(  Vector3::new(0.0, 1.0,  0.0), Vector3::new(-params[0], 0.0, 0.0)),
                        Screw::from_angular_linear(  Vector3::new(0.0, 1.0,  0.0), Vector3::new(-params[0], 0.0, params[2])),
                        Screw::from_angular_linear(  Vector3::new(0.0, 1.0,  0.0), Vector3::new(-params[0], 0.0, params[2] + params[3])),
                        Screw::from_angular_linear(  Vector3::new(0.0, 0.0, -1.0), Vector3::new(-params[4], params[2] + params[3], 0.0)),
                        Screw::from_angular_linear(  Vector3::new(0.0, 1.0,  0.0), Vector3::new(params[1] - params[0], 0.0, params[2] + params[3]))
                    ],
                end_effector_at_initial_position: Matrix4::<f64>::from_row_slice(&[
                       -1.0, 0.0, 0.0, params[2] + params[3],
                        0.0, 0.0, 1.0, params[4] + params[5],
                        0.0, 1.0, 0.0, params[0] - params[1],
                        0.0, 0.0, 0.0,     1.0,
                ])
            })
    }

}

impl ScrewChainLike for Revolute3PrismaticRevolute3 {

    fn from_default() -> Option<ScrewChain> {
        let default_parameters: HashMap<&str, f64> = [
            ("l1", 1.0),
            ("l2", 1.0),
            ].iter().cloned().collect();
        Self::from_parameters(&default_parameters)
    }

    fn parameters() -> Vec<String> {
        return vec![String::from("l1"), String::from("l2")];
    }

    fn from_parameters(parameters_map: &HashMap<&str, f64>) -> Option<ScrewChain> {
        let params_opt = Self::param_map_to_param_vect(parameters_map);
        params_opt.map(|params|
            ScrewChain {
                screws: vec![
                        Screw::from_angular_linear(  Vector3::new( 0.0, 0.0, 1.0), Vector3::new(0.0, 0.0, 0.0)),
                        Screw::from_angular_linear(  Vector3::new( 1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0)),
                        Screw::from_angular_linear(  Vector3::new( 0.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)),
                        Screw::from_angular_linear(  Vector3::new( 0.0, 1.0, 0.0), Vector3::new(0.0, 0.0, 0.0)),
                        Screw::from_angular_linear(  Vector3::new( 1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, -params[0])),
                        Screw::from_angular_linear(  Vector3::new( 0.0, 1.0, 0.0), Vector3::new(0.0, 0.0, 0.0))
                    ],
                end_effector_at_initial_position: Matrix4::<f64>::from_row_slice(&[
                        1.0, 0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0, params[0] + params[1],
                        0.0, 0.0, 1.0, 0.0,
                        0.0, 0.0, 0.0, 1.0,
                ])
            })
    }

}



impl ScrewChainLike for Revolute3Prismatic {

    fn from_default() -> Option<ScrewChain> {
        let default_parameters: HashMap<&str, f64> = [
            ("l1", 1.0),
            ("l2", 1.0),
            ].iter().cloned().collect();
        Self::from_parameters(&default_parameters)
    }

    fn parameters() -> Vec<String> {
        return vec![String::from("l1"), String::from("l2")];
    }

    fn from_parameters(parameters_map: &HashMap<&str, f64>) -> Option<ScrewChain> {
        let params_opt = Self::param_map_to_param_vect(parameters_map);
        params_opt.map(|params|
            ScrewChain {
                screws: vec![
                        Screw::from_angular_linear(  Vector3::new( 0.0, 0.0, 1.0), Vector3::new(0.0, 0.0,  0.0)),
                        Screw::from_angular_linear(  Vector3::new( 0.0, 0.0, 1.0), Vector3::new(0.0, -params[0], 0.0)),
                        Screw::from_angular_linear(  Vector3::new( 0.0, 0.0, 1.0), Vector3::new(0.0, -params[0]-params[1], 0.0)),
                        Screw::from_angular_linear(  Vector3::new( 0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0)),
                    ],
                end_effector_at_initial_position: Matrix4::<f64>::from_row_slice(&[
                    1.0, 0.0, 0.0, params[0] + params[1],
                    0.0, 1.0, 0.0, 0.0,
                    0.0, 0.0, 1.0, 0.0,
                    0.0, 0.0, 0.0, 1.0,
                ])
            })
    }

}
