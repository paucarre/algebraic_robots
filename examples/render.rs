use bevy::prelude::*;
use bevy::input::mouse::{MouseWheel,MouseMotion};
use bevy::render::camera::PerspectiveProjection;
use bevy_prototype_debug_lines::*;
use nalgebra::{ Vector3, Matrix4};
use algebraic_robots::algebraic_robots::*;
use algebraic_robots::screw_chains::*;

//UniversalRobotsUR5
struct PanOrbitCamera {
    /// Mouse
    pub focus: Vec3,
    pub radius: f32,
    pub upside_down: bool,
    pub configuration: CameraConfiguration,
    pub linear_kinematic_state: KinematicState,
    pub rotational_kinematic_state: KinematicState,
    pub camera_position: Vec3,
    pub camera_orientation: Quat,
}

struct CameraConfiguration {
    // Keyboard
	pub linear_keys: Keys,
    pub rotation_keys: Keys,
	/// If `false`, disable keyboard control of the camera. Defaults to `true`
	pub enabled: bool,
    pub key_rotation: KeyCode,
}

struct KinematicState {
   acceleration: f32,
   max_speed: f32,
   friction: f32,
   velocity: Vec3,
}

struct Keys {
    key_left: KeyCode,
    key_right: KeyCode,
    key_up: KeyCode,
    key_down: KeyCode,
    key_forward: KeyCode,
    key_backward: KeyCode,
}

trait KinematicStateOperations {
    fn compute_acceleration(&self, force: Vec3) -> Vec3;
    fn compute_friction(&self) -> Vec3;
    fn compute_velocity(&self, friction: Vec3, time: &Res<Time>) -> Vec3;
    fn apply_acceleration(&mut self, time: &Res<Time>, force: Vec3);
    fn update_kinematic_state(&mut self, time: &Res<Time>);
    fn as_rotation(&self) -> Quat;
    fn as_translation(&self) -> Vec3;
}

impl KinematicStateOperations for KinematicState {


    fn compute_acceleration(&self, force: Vec3) -> Vec3 {
        let acceleration = force;
        if acceleration.length() != 0.0 {
            acceleration.normalize() * self.acceleration
        } else {
            Vec3::ZERO
        }
    }

    fn compute_friction(&self) -> Vec3 {
        if self.velocity.length() != 0.0 {
            self.velocity.normalize() * -1.0 * self.friction
        } else {
            Vec3::ZERO
        }
    }

    fn compute_velocity(&self, friction: Vec3, time: &Res<Time>) -> Vec3 {
        let velocity_clamped = if self.velocity.length() > self.max_speed {
            self.velocity.normalize() * self.max_speed
        } else {
            self.velocity
        };
        let delta_friction = friction * time.delta_seconds();
        let sign_difference = (velocity_clamped.abs() - delta_friction.abs()).signum() / 2.0;
        (velocity_clamped + delta_friction) * (Vec3::new(0.5, 0.5, 0.5) + sign_difference)
    }

    fn apply_acceleration(&mut self, time: &Res<Time>, force: Vec3) {
        let acceleration = self.compute_acceleration(force);
		self.velocity += acceleration * time.delta_seconds();
    }

    fn update_kinematic_state(&mut self, time: &Res<Time>) {
		let friction: Vec3 = self.compute_friction();
        self.velocity = self.compute_velocity(friction, &time);
    }

    fn as_rotation(&self) -> Quat {
        let pitch = Quat::from_rotation_x(self.velocity.x);
        let yaw = Quat::from_rotation_y(self.velocity.y);
        let roll = Quat::from_rotation_z(self.velocity.z);
        pitch * yaw * roll        
    }

    fn as_translation(&self) -> Vec3 {
        Vec3::new(self.velocity.x, self.velocity.y, self.velocity.z)
    }
}

impl Default for PanOrbitCamera {
    fn default() -> Self {
        PanOrbitCamera {
            // Mouse 3D orbit and pan
            focus: Vec3::ZERO,
            radius: 5.0,
            upside_down: false,
            configuration: CameraConfiguration {
                // Keybaord 3D PAN
                key_rotation: KeyCode::LShift,
                linear_keys: Keys {
                    key_left: KeyCode::A,
                    key_right: KeyCode::D,
                    key_up: KeyCode::R,
                    key_down: KeyCode::F,
                    key_forward: KeyCode::W,
                    key_backward: KeyCode::S,
                },
                rotation_keys: Keys {
                    key_left: KeyCode::R,
                    key_right: KeyCode::F,
                    key_up: KeyCode::A,
                    key_down: KeyCode::D,
                    key_forward: KeyCode::W,
                    key_backward: KeyCode::S,
                },
                enabled: true,
            },
            rotational_kinematic_state: KinematicState {
                acceleration: 0.02,
                max_speed: 0.4,
                friction: 0.01,
                velocity: Vec3::ZERO,
            },
            linear_kinematic_state: KinematicState {
                acceleration: 0.2,
                max_speed: 0.5,
                friction: 0.15,
                velocity: Vec3::ZERO,
            },
            camera_position: Vec3::new(0.0, 0.0, 10.0),
            camera_orientation: Quat::from_xyzw(0.0, 0.0, 0.0, 1.0),
        }
    }
}

pub fn movement_axis(
	input: &Res<Input<KeyCode>>,
	plus: KeyCode,
	minus: KeyCode,
) -> f32 {
	let mut axis = 0.0;
	if input.pressed(plus) {
		axis += 1.0;
	}
	if input.pressed(minus) {
		axis -= 1.0;
	}
	axis
}

fn get_key_force(is_enabled: bool, keyboard_input: &Res<Input<KeyCode>>, keys: &Keys) -> Vec3 {
    if is_enabled {
        Vec3::new(
            movement_axis(&keyboard_input, keys.key_right   , keys.key_left),
            movement_axis(&keyboard_input, keys.key_up      , keys.key_down),
            movement_axis(&keyboard_input, keys.key_backward, keys.key_forward),
        )
    } else {
        Vec3::new(0.0, 0.0, 0.0)
    }
}

fn keyboard_camera_system(
	time: Res<Time>,
	keyboard_input: Res<Input<KeyCode>>,
    mut query: Query<(&mut PanOrbitCamera, &mut Transform)>,
) {
	for (mut camera, mut transform) in query.iter_mut() {
        let is_rotation = keyboard_input.pressed(camera.configuration.key_rotation);
        let keys = if is_rotation {
            &camera.configuration.rotation_keys
        } else {
            &camera.configuration.linear_keys
        };
        let force = get_key_force(camera.configuration.enabled, &keyboard_input, &keys);
        if is_rotation {
            camera.rotational_kinematic_state.apply_acceleration(&time, force);
        } else {
            camera.linear_kinematic_state.apply_acceleration(&time, force);
        }
        camera.rotational_kinematic_state.update_kinematic_state(&time);
        let rotation = camera.rotational_kinematic_state.as_rotation();
        camera.camera_orientation = camera.camera_orientation * rotation; 
        camera.linear_kinematic_state.update_kinematic_state(&time);
        let translation = camera.linear_kinematic_state.as_translation();
        camera.camera_position += translation;
        transform.rotation = camera.camera_orientation;
        let camera_rotation = Mat3::from_quat(camera.camera_orientation);
        transform.translation = camera_rotation.mul_vec3(camera.camera_position);
 	}
}

//fn robot_draw_system(time: Res<Time>, mut lines: ResMut<DebugLines>) {
    //let seconds = time.seconds_since_startup() as f32;
    //lines.line_colored(Vec3::new(0.0, 0.0, 0.0), Vec3::new(f32::sin(seconds + 3.14), 1.0, 0.0),  0.0, Color::WHITE);
//}

pub struct RobotDraw {
    pub screw_chain: ScrewChain,
    pub coordinates : Vec<f64>,
}

impl Default for RobotDraw {
    fn default() -> Self {
        Self {
            screw_chain: UniversalRobotsUR5::from_default().unwrap(),
            coordinates: vec![],
        }
    }
}

impl RobotDraw {
    pub fn apply_coordinates(&mut self, coordinates : &[f64], duration: f32) {
      self.coordinates = coordinates.to_vec();
    }
}


fn robot_draw_system(
    mut lines: ResMut<DebugLines>,
    mut robot: ResMut<RobotDraw>,
    time: Res<Time>,
) {
    let seconds = time.seconds_since_startup() as f32;
    lines.line_colored(Vec3::new(0.0, 0.0, 0.0), Vec3::new(f32::sin(seconds + 3.14), 1.0, 0.0),  0.0, Color::WHITE);
}



/// Pan the camera with middle mouse click, zoom with scroll wheel, orbit with right mouse click.
fn mouse_camera_system(
    windows: Res<Windows>,
    mut ev_motion: EventReader<MouseMotion>,
    mut ev_scroll: EventReader<MouseWheel>,
    input_mouse: Res<Input<MouseButton>>,
    mut query: Query<(&mut PanOrbitCamera, &mut Transform, &PerspectiveProjection)>,
) {
    // change input mapping for orbit and panning here
    let orbit_button = MouseButton::Right;
    let pan_button = MouseButton::Middle;

    let mut pan = Vec2::ZERO;
    let mut rotation_move = Vec2::ZERO;
    let mut scroll = 0.0;
    let mut orbit_button_changed = false;

    if input_mouse.pressed(orbit_button) {
        for ev in ev_motion.iter() {
            rotation_move += ev.delta;
        }
    } else if input_mouse.pressed(pan_button) {
        // Pan only if we're not rotating at the moment
        for ev in ev_motion.iter() {
            pan += ev.delta;
        }
    }
    for ev in ev_scroll.iter() {
        scroll += ev.y;
    }
    if input_mouse.just_released(orbit_button) || input_mouse.just_pressed(orbit_button) {
        orbit_button_changed = true;
    }

    for (mut pan_orbit, mut transform, projection) in query.iter_mut() {
        if orbit_button_changed {
            // only check for upside down when orbiting started or ended this frame
            // if the camera is "upside" down, panning horizontally would be inverted, so invert the input to make it correct
            let up = transform.rotation * Vec3::Y;
            pan_orbit.upside_down = up.y <= 0.0;
        }

        let camera_is_updated = if rotation_move.length_squared() > 0.0 {
            let window = get_primary_window_size(&windows);
            let delta_x = {
                let delta = rotation_move.x / window.x * std::f32::consts::PI * 2.0;
               if pan_orbit.upside_down { -delta } else { delta }
            };
            let delta_y = rotation_move.y / window.y * std::f32::consts::PI;
            let yaw = Quat::from_rotation_y(-delta_x);
            let pitch = Quat::from_rotation_x(-delta_y);
            transform.rotation = yaw * transform.rotation; // rotate around global y axis
            transform.rotation = transform.rotation * pitch; // rotate around local x axis
            true
        } else if pan.length_squared() > 0.0 {
            // make panning distance independent of resolution and FOV,
            let window = get_primary_window_size(&windows);
            pan *= Vec2::new(projection.fov * projection.aspect_ratio, projection.fov) / window;
            // translate by local axes
            let right = transform.rotation * Vec3::X * -pan.x;
            let up = transform.rotation * Vec3::Y * pan.y;
            // make panning proportional to distance away from focus point
            let translation = (right + up) * pan_orbit.radius;
            pan_orbit.focus += translation;
            true
        } else if scroll.abs() > 0.0 {
            pan_orbit.camera_position.z -= scroll * pan_orbit.radius * 0.2;
            // dont allow zoom to reach zero or you get stuck
            pan_orbit.radius = f32::max(pan_orbit.radius, 0.05);
            true
        } else {
            false
        };

        if camera_is_updated {
            // emulating parent/child to make the yaw/y-axis rotation behave like a turntable
            // parent = x and y rotation
            // child = z-offset
            let rot_matrix = Mat3::from_quat(transform.rotation);
            transform.translation = rot_matrix.mul_vec3(pan_orbit.camera_position);
        }
    }
}

fn get_primary_window_size(windows: &Res<Windows>) -> Vec2 {
    let window = windows.get_primary().unwrap();
    let window = Vec2::new(window.width() as f32, window.height() as f32);
    window
}

/// Spawn a camera like this
fn spawn_camera(mut commands: Commands) {
    let translation = Vec3::new(-2.0, 2.5, 5.0);
    let radius = translation.length();

    commands.spawn_bundle(PerspectiveCameraBundle {
        transform: Transform::from_translation(translation)
            .looking_at(Vec3::ZERO, Vec3::Y),
        ..Default::default()
    }).insert(PanOrbitCamera {
        radius,
        ..Default::default()
    });
}




// ANCHOR_END: example

fn spawn_scene(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    // spawn a cube and a light
    commands.spawn_bundle(PbrBundle {
        mesh: meshes.add(Mesh::from(shape::Cube { size: 1.0 })),
        material: materials.add(Color::rgb(0.8, 0.7, 0.6).into()),
        transform: Transform::from_translation(Vec3::new(0.0, 0.5, 0.0)),
        ..Default::default()
    });
    commands.spawn_bundle(LightBundle {
        transform: Transform::from_translation(Vec3::new(4.0, 8.0, 4.0)),
        ..Default::default()
    });

    spawn_camera(commands);
}

fn main() {
    App::build()
        .init_resource::<RobotDraw>()
        .add_plugins(DefaultPlugins)
        .add_plugin(DebugLinesPlugin)
        .add_startup_system(spawn_scene.system())
        .add_system(mouse_camera_system.system())
        .add_system(keyboard_camera_system.system())
        .add_system_to_stage(CoreStage::Last, robot_draw_system.system().label("robot_draw"))
        .run();

    // just to catch compilation errors
    let _ = App::build()
        .add_startup_system(spawn_camera.system());
}