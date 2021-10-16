use bevy::prelude::*;
use bevy::input::mouse::{MouseWheel,MouseMotion};
use bevy::render::camera::PerspectiveProjection;

// ANCHOR: example
/// Tags an entity as capable of panning and orbiting.
struct PanOrbitCamera {
    /// The "focus point" to orbit around. It is automatically updated when panning the camera
    pub focus: Vec3,
    pub radius: f32,
    pub upside_down: bool,
    //2D
    // The speed the FlyCamera2d accelerates at.
	/// The maximum speed the FlyCamera can move at.
	/// The amount of deceleration to apply to the camera's motion.
	/// The current velocity of the FlyCamera2d. This value is always up-to-date, enforced by [FlyCameraPlugin](struct.FlyCameraPlugin.html)
	/// Key used to move left. Defaults to <kbd>A</kbd>
	pub key_left: KeyCode,
	/// Key used to move right. Defaults to <kbd>D</kbd>
	pub key_right: KeyCode,
	/// Key used to move up. Defaults to <kbd>W</kbd>
	pub key_up: KeyCode,
	/// Key used to move forward. Defaults to <kbd>S</kbd>
	pub key_down: KeyCode,
	/// If `false`, disable keyboard control of the camera. Defaults to `true`
	pub enabled: bool,
    pub key_forward: KeyCode,
    pub key_backward: KeyCode,
    pub key_rotation: KeyCode,
    pub linear_kinematic_state: KinematicState,
    pub rotational_kinematic_state: KinematicState,
}

struct KinematicState {
   acceleration: f32,
   max_speed: f32,
   friction: f32,
   velocity: Vec3,
}

trait KinematicStateOperations {
    fn compute_acceleration(&self, 
        axis_horizontal: f32, axis_vertical: f32, axis_depth: f32) -> Vec3;
    fn compute_friction(&self) -> Vec3;
    fn compute_velocity(&self, friction: Vec3, time: &Res<Time>) -> Vec3;
    fn apply_acceleration(&mut self, time: &Res<Time>, axis_horizontal: f32, axis_vertical: f32, axis_depth: f32);
    fn update_kinematic_state(&mut self, time: &Res<Time>);
    fn as_rotation(&self) -> Quat;
    fn as_translation(&self) -> Vec3;
}

impl KinematicStateOperations for KinematicState {


    fn compute_acceleration(&self, 
        axis_horizontal: f32, axis_vertical: f32, axis_depth: f32) -> Vec3 {
        let acceleration: Vec3 = (Vec3::X * axis_horizontal) + 
            (Vec3::Y * axis_vertical) + 
            (Vec3::Z * axis_depth);
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

    fn apply_acceleration(&mut self, time: &Res<Time>, axis_horizontal: f32, axis_vertical: f32, axis_depth: f32) {
        let acceleration = self.compute_acceleration(axis_horizontal, axis_vertical, axis_depth);
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
		const MUL_2D: f32 = 10.0;
        PanOrbitCamera {
            // Mouse 3D orbit and pan
            focus: Vec3::ZERO,
            radius: 5.0,
            upside_down: false,
            // Keybaord 3D PAN
            key_left: KeyCode::A,
            key_right: KeyCode::D,
            key_up: KeyCode::R,
            key_down: KeyCode::F,
            key_forward: KeyCode::W,
            key_backward: KeyCode::S,
            key_rotation: KeyCode::LShift,
            enabled: true,
            rotational_kinematic_state: KinematicState {
                acceleration: 0.02,
                max_speed: 0.4,
                friction: 0.01,
                velocity: Vec3::ZERO,
            },
            linear_kinematic_state: KinematicState {
                acceleration: 0.2,
                max_speed: 1.0,
                friction: 0.1,
                velocity: Vec3::ZERO,
            }
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


fn camera_2d_movement_system(
	time: Res<Time>,
	keyboard_input: Res<Input<KeyCode>>,
    mut query: Query<(&mut PanOrbitCamera, &mut Transform)>,
) {
	for (mut options, mut transform) in query.iter_mut() {
		let (axis_horizontal, axis_vertical, axis_depth) = if options.enabled {
			(
				movement_axis(&keyboard_input, options.key_right, options.key_left),
				movement_axis(&keyboard_input, options.key_up, options.key_down),
				movement_axis(&keyboard_input, options.key_backward, options.key_forward),
			)
		} else {
			(0.0, 0.0, 0.0)
		};

        if keyboard_input.pressed(options.key_rotation) {
            options.rotational_kinematic_state.apply_acceleration(&time, axis_horizontal, axis_vertical, axis_depth);
        } else {
            options.linear_kinematic_state.apply_acceleration(&time, axis_horizontal, axis_vertical, axis_depth);
        }
        options.rotational_kinematic_state.update_kinematic_state(&time);
        let rotation = options.rotational_kinematic_state.as_rotation();
        transform.rotation = transform.rotation * rotation;
        options.linear_kinematic_state.update_kinematic_state(&time);
        let translation = options.linear_kinematic_state.as_translation();
        let rot_matrix = Mat3::from_quat(transform.rotation);
        transform.translation += rot_matrix.mul_vec3(translation);
	}
}

/// Pan the camera with middle mouse click, zoom with scroll wheel, orbit with right mouse click.
fn pan_orbit_camera(
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

        let mut any = false;
        if rotation_move.length_squared() > 0.0 {
            any = true;
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
        } else if pan.length_squared() > 0.0 {
            any = true;
            // make panning distance independent of resolution and FOV,
            let window = get_primary_window_size(&windows);
            pan *= Vec2::new(projection.fov * projection.aspect_ratio, projection.fov) / window;
            // translate by local axes
            let right = transform.rotation * Vec3::X * -pan.x;
            let up = transform.rotation * Vec3::Y * pan.y;
            // make panning proportional to distance away from focus point
            let translation = (right + up) * pan_orbit.radius;
            pan_orbit.focus += translation;
        } else if scroll.abs() > 0.0 {
            any = true;
            pan_orbit.radius -= scroll * pan_orbit.radius * 0.2;
            // dont allow zoom to reach zero or you get stuck
            pan_orbit.radius = f32::max(pan_orbit.radius, 0.05);
        }

        if any {
            // emulating parent/child to make the yaw/y-axis rotation behave like a turntable
            // parent = x and y rotation
            // child = z-offset
            let rot_matrix = Mat3::from_quat(transform.rotation);
            transform.translation = rot_matrix.mul_vec3(Vec3::new(0.0, 0.0, pan_orbit.radius));
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
        .add_plugins(DefaultPlugins)
        .add_startup_system(spawn_scene.system())
        .add_system(pan_orbit_camera.system())
        .add_system(camera_2d_movement_system.system())
        .run();

    // just to catch compilation errors
    let _ = App::build()
        .add_startup_system(spawn_camera.system());
}