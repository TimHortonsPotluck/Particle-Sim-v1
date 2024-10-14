#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <SFML/Graphics.hpp>

// hue: 0-360°; sat: 0.f-1.f; val: 0.f-1.f
const sf::Color hsv(int hue, float sat, float val) {
	hue %= 360;
	while (hue < 0) hue += 360;

	if (sat < 0.f) sat = 0.f;
	if (sat > 1.f) sat = 1.f;

	if (val < 0.f) val = 0.f;
	if (val > 1.f) val = 1.f;

	int h = hue / 60;
	float f = float(hue) / 60 - h;
	float p = val * (1.f - sat);
	float q = val * (1.f - sat * f);
	float t = val * (1.f - sat * (1 - f));

	switch (h) {
	default:
	case 0:
	case 6: return sf::Color(val * 255, t * 255, p * 255);
	case 1: return sf::Color(q * 255, val * 255, p * 255);
	case 2: return sf::Color(p * 255, val * 255, t * 255);
	case 3: return sf::Color(p * 255, q * 255, val * 255);
	case 4: return sf::Color(t * 255, p * 255, val * 255);
	case 5: return sf::Color(val * 255, p * 255, q * 255);
	}
}

class BeemanObject {

public:
	sf::Vector2<double> pos;
	sf::Vector2<double> next_pos;
	sf::Vector2<double> vel;
	sf::Vector2<double> next_vel_predicted;
	sf::Vector2<double> acc;
	sf::Vector2<double> last_acc;
	sf::Vector2<double> next_acc;
	double radius = 10;
	sf::Color color = sf::Color::Cyan;
	unsigned int id = -1;
	unsigned int counter = 0;
	double pressure = 0;

	BeemanObject(sf::Vector2<double> _pos, sf::Vector2<double> _vel, double _radius) : 
		pos{ _pos },
		next_pos{ _pos },
		vel{ _vel },
		next_vel_predicted{ _vel },
		acc{ sf::Vector2<double>(0, 0) }, 
		last_acc{ sf::Vector2<double>(0, 0) },
		next_acc{ sf::Vector2<double>(0, 0) },
		radius{ _radius } {}

	void predictNextState(double dt) {
		next_pos = pos + vel * dt + (4. * acc - last_acc) * dt * dt / 6.;
		next_vel_predicted = vel + (3. * acc - last_acc) * .5 * dt; // predicted

	}

	void update(double dt) {

		pos = next_pos;
		const sf::Vector2<double> next_vel_corrected = vel + (5. * next_acc + 8. * acc - last_acc) * dt / 12.; // the "real" vel
		vel = next_vel_corrected;
		last_acc = acc;
		acc = next_acc;
		next_acc = { 0, 0 };

	}

	void addNextAcc(sf::Vector2<double> a) {
		next_acc += a;
	}

	void setVel(sf::Vector2<double> v) {
		vel = v;
	}

	double getVelMagSq() {
		return vel.x * vel.x + vel.y * vel.y;
	}

	double getVelMag() {
		return sqrt(vel.x * vel.x + vel.y * vel.y);
	}
};

class PhysicsObject {

public:
	const double pi = 3.14159265358979323846;
	const double to_degrees = 180 / pi;
	sf::Vector2<double> pos;
	double angle;
	sf::Vector2<double> last_pos;
	double last_angle;
	sf::Vector2<double> acc;
	double ang_acc;
	double radius = 10;
	sf::Color color = sf::Color::White;
	unsigned int id = -1;
	unsigned int counter = 0;
	double forces = 0;
	double forces_color_amt = 1;
	double charge;

	PhysicsObject(sf::Vector2<double> _pos, double _angle, double _charge, double _radius) : pos{ _pos }, angle{ _angle }, last_pos{ _pos }, last_angle{ _angle }, acc{ sf::Vector2<double>(0, 0) }, charge{ _charge }, ang_acc { 0 }, radius{ _radius } {}

	PhysicsObject(sf::Vector2<double> _pos, double _angle, double _radius) : pos{ _pos }, angle{ _angle }, last_pos{ _pos }, last_angle{ _angle }, acc{ sf::Vector2<double>(0, 0) }, charge{ 0 }, ang_acc{ 0 }, radius{ _radius } {}

	void update(double dt) {
		const sf::Vector2<double> next_pos = 2. * pos - last_pos + acc * dt * dt;
		last_pos = pos;
		pos = next_pos;
		acc = { 0, 0 };

		const double next_angle = 2. * angle - last_angle + ang_acc * dt * dt;
		last_angle = angle;
		angle = next_angle;
		ang_acc = 0;

		//color = sf::Color::Color(255, std::max(0., 255 - forces / forces_color_amt), std::max(0., 255 - forces / forces_color_amt)); // the coloring step should be done via renderer
		forces = 0;
	}

	void addAcc(sf::Vector2<double> a) {
		acc += a;
		forces += sqrt(a.x * a.x + a.y * a.y);
	}

	void setVel(sf::Vector2<double> v, double dt) {
		last_pos = pos - v * dt;
	}
	
	void addVel(sf::Vector2<double> v, double dt) {
		last_pos -= v * dt;
	}

	sf::Vector2<double> getVel(double dt) const {
		return (pos - last_pos) / dt;
	}

	void addAngAcc(double a_a) {
		ang_acc += a_a;
	}

	void setAngVel(double a_v, double dt) {
		last_angle = angle - a_v * dt;
	}

	void addAngVel(double a_v, double dt) {
		last_angle -= a_v * dt;
	}

	double getAngVel(double dt) const {
		return (angle - last_angle) / dt;
	}

	void resetForcesPressure() {
		forces = 0;
	}

	const double getVelMagSq(double dt) const {
		sf::Vector2<double> vel = getVel(dt);
		return vel.x * vel.x + vel.y * vel.y;
	}

	const double getVelMag(double dt) const {
		sf::Vector2<double> vel = getVel(dt);
		return sqrt(vel.x * vel.x + vel.y * vel.y);
	}

	const sf::Color getColorDipole() const {
		return hsv(angle * to_degrees, 1, 1);
	}

	const sf::Color getColorCharge() const {
		return sf::Color::Color(std::max(charge * 255, 0.), std::max(-charge * 255, 0.), 0);
	}
};

class Solver {

public:
	Solver(sf::Vector2u _window_size, bool _do_dipole) : window_size{ _window_size }, do_dipole{ _do_dipole } {};

	void addObject(PhysicsObject& obj) {
		obj.id = getObjectsCount();
		objects.push_back(obj);
	}

	void addObject(PhysicsObject obj) {
		obj.id = getObjectsCount();
		objects.push_back(obj);
	}

	void addObject(sf::Vector2<double> pos, double angle, double radius) {
		//obj.id = getObjectsCount();
		objects.emplace_back(pos, angle, radius);
	}

	void addObject(sf::Vector2<double> pos, double angle, double charge, double radius) {
		//obj.id = getObjectsCount();
		objects.emplace_back(pos, angle, charge, radius);
	}

	void update() {
		time += frame_dt;
		const double step_dt = getStepDT();
		//std::cout << step_dt << std::endl;
		for (int i = 0; i < sub_steps; i++) {
			//addGravity();
			//addNewLennardJones();
			//addDipoleDipoleForce();
			addWallRepellant();
			if (do_dipole) {
				addCombinedLJDipoleForce();
			} else {
				addNewLennardJones();
			}
			//addCombinedLJChargeForce();
			//addAngleAirFriction(step_dt);
			addAirFriction(step_dt);
			//addDebugForce();
			//addCenterGravity();
			updateObjects(step_dt);
		}
	}

	void setSimUpdateRate(int rate) {
		frame_dt = 1 / static_cast<double>(rate);
		//std::cout << frame_dt << std::endl;
	}

	void setSubStepsCount(int _sub_steps) {
		sub_steps = _sub_steps;
	}

	void setObjVel(PhysicsObject& obj, sf::Vector2<double> v) {
		obj.setVel(v, getStepDT());
	}

	const std::vector<PhysicsObject>& getObjects() const {
		return objects;
	}

	size_t getObjectsCount() const {
		return objects.size();
	}

	double getTime() const {
		return time;
	}

	double getStepDT() const {
		return frame_dt / static_cast<double>(sub_steps);
	}

	void setGravity(sf::Vector2<double> g) {
		gravity = g;
	}

	double getTotalVelSq(double dt) {
		double total = 0;
		for (PhysicsObject& obj : objects) {
			total += obj.getVelMagSq(dt);
		}
		return total;
	}

	double calcTotalEnergyLJ_debug() { // 2-particle system assumed
		
		sf::Vector2<double> disp = objects[1].pos - objects[0].pos;
		const double dist_sq = disp.x * disp.x + disp.y * disp.y;
		const double inv_dist_sq = 1. / dist_sq;
		const double b = (lj_strength / 6.) * (lj_strength * inv_dist_sq) * (lj_strength * inv_dist_sq) * (lj_strength * inv_dist_sq);
		const double z = lj_r0 * lj_r0 / dist_sq;
		const sf::Vector2<double> v = objects[0].getVel(getStepDT());

		//std::cout << disp.x << std::endl;
		//std::cout << (lj_strength / 6.) << std::endl;
		//std::cout << (lj_strength * inv_dist_sq) << std::endl;
		//std::cout << (lj_strength * inv_dist_sq)* (lj_strength * inv_dist_sq)* (lj_strength * inv_dist_sq) << std::endl;
		//std::cout << (lj_strength * inv_dist_sq) * (lj_strength * inv_dist_sq) * (lj_strength * inv_dist_sq)* lj_strength/6 << std::endl;
		//std::cout << b << std::endl;
		const double vmag = v.x * v.x + v.y * v.y;
		const double amt = (b * (.5 * z * z * z - 1));
		//std::cout << amt  << "=" << std::endl;
		return amt + vmag * .5;
	}

	double calcTotalEnergyDebugForce_debug() { // 2-particle system assumed

		sf::Vector2<double> disp = objects[1].pos - objects[0].pos;
		const double dist_sq = disp.x * disp.x + disp.y * disp.y;
		const double dist = sqrt(dist_sq);
		const sf::Vector2<double> v = objects[0].getVel(getStepDT());
		//std::cout << v.x << std::endl;
		const double vmag = v.x * v.x + v.y * v.y;
		return debug_force_strength * ((1 / dist_sq) - (1 / dist));// +vmag * .5;
	}

	double calcTotalEnergyGravity_debug() { // only for single particle

		const double y = 400 - objects[0].pos.y;
		const sf::Vector2<double> v = objects[0].getVel(getStepDT());
		//std::cout << v.x << std::endl;
		const double vmag = v.x * v.x + v.y * v.y;
		return y * gravity.y + vmag * .5;
	}

	double calcTotalEnergyCenterGravity_debug() {
		sf::Vector2<double> center = { 200, 200 };
		PhysicsObject& this_obj = objects[0];
		const sf::Vector2<double> disp = this_obj.pos - center;
		const double dist_sq = disp.x * disp.x + disp.y * disp.y;
		const double dist = sqrt(dist_sq);
		const sf::Vector2<double> v = objects[0].getVel(getStepDT());
		const double vmag = v.x * v.x + v.y * v.y;
		return -center_grav_strength / dist;// +vmag * .5;
	}

	void setLJR0(double _lj_r0) {
		lj_r0 = _lj_r0;
		lj_A = lj_epsilon * pow(lj_r0, 12);
		lj_B = 2 * lj_epsilon * pow(lj_r0, 6);
		lj_A_wall = lj_epsilon * pow(lj_r0 * .5, 12);
	}

	void setLJEpsilon(double _lj_epsilon) {
		lj_epsilon = _lj_epsilon;
		lj_A = lj_epsilon * pow(lj_r0, 12);
		lj_B = 2 * lj_epsilon * pow(lj_r0, 6);
		lj_A_wall = lj_epsilon * pow(lj_r0 * .5, 12);
	}

	void setDipoleK(double _dipole_k) {
		dipole_k = _dipole_k;
	}

	void setAirFrictionCoeff(double _air_friction_coeff) {
		air_friction_coeff = _air_friction_coeff;
	}

	void setAirAngleFrictionCoeff(double _air_angle_friction_coeff) {
		air_angle_friction_coeff = _air_angle_friction_coeff;
	}

	void setDoDipole(bool _do_dipole) {
		do_dipole = _do_dipole;
	}

	bool getDoDipole() {
		return do_dipole;
	}

private:
	int sub_steps = 1;
	sf::Vector2<double> gravity = {0, .1};

	sf::Vector2u window_size;
	bool do_dipole;

	const double lj_strength = 500; // The "actual" factor that wil multiply the LJ force will be this value to the 4th power, for computation reasons. Not used in new LJ

	double lj_r0 = 20; // The equilibrium radius for the LJ potential.

	double lj_epsilon = 10; // for the "new" LJ 
	double lj_A = lj_epsilon * pow(lj_r0, 12);
	double lj_B = 2 * lj_epsilon * pow(lj_r0, 6);
	double lj_A_wall = lj_epsilon * pow(lj_r0 *.5, 12);
	
	double dipole_k = 100000;

	const double charge_k = 1000;

	const double debug_force_strength = 10;
	const double center_grav_strength = 1000;
	double air_friction_coeff = .05;
	double air_angle_friction_coeff = .2;
	
	std::vector<PhysicsObject> objects;
	double time = 0;
	double frame_dt = 0;

	void calcNextForces(double dt) {

	}

	void addGravity() {
		for (PhysicsObject& obj : objects) {
			obj.addAcc(gravity);
		}
	}

	void addAirFriction(double dt) {
		for (PhysicsObject& obj : objects) {
			obj.addAcc(-air_friction_coeff * obj.getVel(dt));
		}
	}

	void addAngleAirFriction(double dt) {
		for (PhysicsObject& obj : objects) {
			obj.addAngAcc(-air_angle_friction_coeff * obj.getAngVel(dt));
		}
	}

	void addWallRepellant() {
		for (PhysicsObject& obj : objects) {
			//left wall
			const double left_dist_13 = obj.pos.x * obj.pos.x * obj.pos.x * obj.pos.x * obj.pos.x * obj.pos.x * obj.pos.x * obj.pos.x * obj.pos.x * obj.pos.x * obj.pos.x * obj.pos.x * obj.pos.x;
			const sf::Vector2<double> f_left = { 12 * lj_A_wall / left_dist_13, 0 };
			//right wall
			const double right_dist = window_size.x - obj.pos.x;
			const double right_dist_13 = right_dist * right_dist * right_dist * right_dist * right_dist * right_dist * right_dist * right_dist * right_dist * right_dist * right_dist * right_dist * right_dist;
			const sf::Vector2<double> f_right = { -12 * lj_A_wall / right_dist_13, 0 };
			
			//top wall
			const double top_dist_13 = obj.pos.y * obj.pos.y * obj.pos.y * obj.pos.y * obj.pos.y * obj.pos.y * obj.pos.y * obj.pos.y * obj.pos.y * obj.pos.y * obj.pos.y * obj.pos.y * obj.pos.y;
			const sf::Vector2<double> f_top = { 0, 12 * lj_A_wall / top_dist_13 };
			//bottom wall
			const double bottom_dist = window_size.y - obj.pos.y;
			const double bottom_dist_13 = bottom_dist * bottom_dist * bottom_dist * bottom_dist * bottom_dist * bottom_dist * bottom_dist * bottom_dist * bottom_dist * bottom_dist * bottom_dist * bottom_dist * bottom_dist;
			const sf::Vector2<double> f_bottom = { 0, -12 * lj_A_wall / bottom_dist_13 };
			const sf::Vector2<double> net_force = f_left + f_right + f_top + f_bottom;
			//std::cout << net_force.x << ", " << net_force.y << std::endl;
			obj.addAcc(net_force);
		}
	}

	void addCenterGravity() {
		int object_count = getObjectsCount();
		sf::Vector2<double> center = { 200, 200 };
		for (int i = 0; i < object_count; i++) {
			PhysicsObject& this_obj = objects[i];
			const sf::Vector2<double> disp = this_obj.pos - center;
			const double dist_sq = disp.x * disp.x + disp.y * disp.y;
			const double dist = sqrt(dist_sq);
			const sf::Vector2<double> force = (1 / dist_sq) * center_grav_strength * (-disp / dist);
			this_obj.addAcc(force);
		}
	}

	void addDebugForce() {
		// potential is 1/x^2 - 1/x, x is displacement
		// force is 2/x^3 - 1/x^2
		int object_count = getObjectsCount();


		for (int i = 0; i < object_count; i++) {
			PhysicsObject& this_obj = objects[i];
			for (int j = i + 1; j < object_count; j++) {
				PhysicsObject& other_obj = objects[j];
				const sf::Vector2<double> disp = this_obj.pos - other_obj.pos;
				const double dist_sq = disp.x * disp.x + disp.y * disp.y;
				const double dist = sqrt(dist_sq);

				const sf::Vector2<double> force = debug_force_strength * (disp / dist) * ((2 / (dist * dist * dist)) - (1 / dist_sq));
				//std::cout << force.x << std::endl;
				this_obj.addAcc(force); // assuming the mass is 1 for all objects for now
				other_obj.addAcc(-force);

			}
		}
	}

	void addNewLennardJones() {
		//
		//	12A		  6B
		//	---	  -   --
		//	r^13	  r^7
		//

		int object_count = getObjectsCount();
		for (int i = 0; i < object_count; i++) {
			PhysicsObject& this_obj = objects[i];
			for (int j = i + 1; j < object_count; j++) {
				PhysicsObject& other_obj = objects[j];
				const sf::Vector2<double> disp = this_obj.pos - other_obj.pos;
				const double dist_sq = disp.x * disp.x + disp.y * disp.y;

				const double inv_dist_sq = 1 / dist_sq;
				const double inv_dist_6th = inv_dist_sq * inv_dist_sq * inv_dist_sq;

				const sf::Vector2<double> force = disp * (12 * lj_A * inv_dist_6th * inv_dist_6th * inv_dist_sq - 6 * lj_B * inv_dist_6th * inv_dist_sq);
				this_obj.addAcc(force); // assuming the mass is 1 for all objects for now
				other_obj.addAcc(-force);

			}
		}

	}

	void addLennardJones() {
		//	LJ force:
		// 
		//	6b*(disp)	 /	 r0^6		\
		//	--------- * |	------  - 1  |
		//	 disp^8		 \	disp^6		/
		//
		// 6b corresponds to the strength of the force.
		// For the displacements to the 6th and 8th powers, the dist_sq will be used, so as to avoid doing a square root.
		// The variable lj_strength is used to signify the 4th root of the strength (6b), which is divided by the disp_sq. The resultant value will be called d.
		// The variable d is raised to the 4th power rather than disp_sq itself, to make the values being multiplied much smaller and less prone to truncation errors.
		// The variable z is equal to ro^2 / dist_sq, and will be cubed. Again, this is to avoid multiplying and dividing massive numbers.
		// The final equation:
		//			
		//	disp * d^4 * (z^3 - 1)
		// 
		// Where d = ((6b)^(1/8) / disp)^8 = (lj_strength / dist_sq)^4
		// and z = (r0 / disp)^2 = r0^2 / dist_sq
		// 
		// For now, it's assumed all objects have the same potentials and properties.
		// 

		int object_count = getObjectsCount();

		for (int i = 0; i < object_count; i++) {
			PhysicsObject& this_obj = objects[i];
			for (int j = i + 1; j < object_count; j++) {
				PhysicsObject& other_obj = objects[j];
				const sf::Vector2<double> disp = this_obj.pos - other_obj.pos;
				const double dist_sq = disp.x * disp.x + disp.y * disp.y;

				const double d = lj_strength / dist_sq;
				const double z = lj_r0 * lj_r0 / dist_sq;
				const double mag = d * d * d * d * (z * z * z - 1);
				this_obj.addAcc(disp * mag); // assuming the mass is 1 for all objects for now
				other_obj.addAcc(-disp * mag);

			}
		}
	}

	void addCombinedLJDipoleForce() {
		int object_count = getObjectsCount();
		for (int i = 0; i < object_count; i++) {
			PhysicsObject& this_obj = objects[i];
			const sf::Vector2<double> this_dir = { sin(this_obj.angle), cos(this_obj.angle) };
			for (int j = i + 1; j < object_count; j++) {
				PhysicsObject& other_obj = objects[j];
				const sf::Vector2<double> disp = this_obj.pos - other_obj.pos;
				const double dist_sq = disp.x * disp.x + disp.y * disp.y;
				const double dist = sqrt(dist_sq);
				const double inv_dist = 1 / dist;
				const double inv_dist_sq = 1 / dist_sq;
				const double inv_dist_6th = inv_dist_sq * inv_dist_sq * inv_dist_sq;

				//const sf::Vector2<double> lj_force = disp * (12 * lj_A * inv_dist_6th * inv_dist_6th * inv_dist_sq - 6 * lj_B * inv_dist_6th * inv_dist_sq);
				//this_obj.addAcc(lj_force); // assuming the mass is 1 for all objects for now
				//other_obj.addAcc(-lj_force);

				const sf::Vector2<double> disp_norm = disp * inv_dist;
				const sf::Vector2<double> theta_dir = { disp_norm.y, -disp_norm.x };

				const sf::Vector2<double> other_dir = { sin(other_obj.angle), cos(other_obj.angle) };

				const double this_dot_r = disp_norm.x * this_dir.x + disp_norm.y * this_dir.y;
				const double other_dot_r = disp_norm.x * other_dir.x + disp_norm.y * other_dir.y;
				const double this_dot_other = this_dir.x * other_dir.x + this_dir.y * other_dir.y;

				const double this_cross_other = this_dir.x * other_dir.y - this_dir.y * other_dir.x;
				const double this_cross_r = this_dir.x * disp_norm.y - this_dir.y * disp_norm.x;
				const double other_cross_r = other_dir.x * disp_norm.y - other_dir.y * disp_norm.x;

				const double r_component = this_dot_other - 3 * this_dot_r * other_dot_r;
				const double theta_component = this_dot_r * other_cross_r + other_dot_r * this_cross_r; // these aren't actually cross products, but since they're already computed it's simpler to do it this way.

				const double dipole_force_factor = dipole_k * inv_dist_sq * inv_dist_sq;
				const double dipole_torque_factor = -dipole_k * .5 * inv_dist_sq * inv_dist;

				const sf::Vector2<double> total_force = dipole_force_factor * (r_component * disp_norm + theta_component * theta_dir) + disp * (12 * lj_A * inv_dist_6th * inv_dist_6th * inv_dist_sq - 6 * lj_B * inv_dist_6th * inv_dist_sq);

				const double torque_on_other = dipole_torque_factor * (3 * this_dot_r * other_cross_r + this_cross_other);
				const double torque_on_this = dipole_torque_factor * (3 * other_dot_r * this_cross_r - this_cross_other);

				this_obj.addAcc(total_force);
				other_obj.addAcc(-total_force); // assuming the dipole moment is 1 for all objects for now
				
				this_obj.addAngAcc(torque_on_this);
				other_obj.addAngAcc(torque_on_other);
			}
		}
	}

	void addCombinedLJChargeForce() {
		int object_count = getObjectsCount();
		for (int i = 0; i < object_count; i++) {
			PhysicsObject& this_obj = objects[i];
			const sf::Vector2<double> this_dir = { sin(this_obj.angle), cos(this_obj.angle) };
			for (int j = i + 1; j < object_count; j++) {
				PhysicsObject& other_obj = objects[j];
				const sf::Vector2<double> disp = this_obj.pos - other_obj.pos;
				const double dist_sq = disp.x * disp.x + disp.y * disp.y;
				const double dist = sqrt(dist_sq);
				const double inv_dist = 1 / dist;
				const double inv_dist_sq = 1 / dist_sq;
				const double inv_dist_6th = inv_dist_sq * inv_dist_sq * inv_dist_sq;

				//const sf::Vector2<double> lj_force = disp * (12 * lj_A * inv_dist_6th * inv_dist_6th * inv_dist_sq - 6 * lj_B * inv_dist_6th * inv_dist_sq);
				//this_obj.addAcc(lj_force); // assuming the mass is 1 for all objects for now
				//other_obj.addAcc(-lj_force);


				const sf::Vector2<double> total_force = disp * (charge_k * this_obj.charge * other_obj.charge * inv_dist_sq * inv_dist + (12 * lj_A * inv_dist_6th * inv_dist_6th * inv_dist_sq - 6 * lj_B * inv_dist_6th * inv_dist_sq));

				this_obj.addAcc(total_force);
				other_obj.addAcc(-total_force); // assuming the mass is 1 for all objects for now
			}
		}
	}

	void addDipoleDipoleForce() {

		int object_count = getObjectsCount();

		for (int i = 0; i < object_count; i++) {
			PhysicsObject& this_obj = objects[i];
			const sf::Vector2<double> this_dir = { sin(this_obj.angle), cos(this_obj.angle) };
			for (int j = i + 1; j < object_count; j++) {
				PhysicsObject& other_obj = objects[j];
				const sf::Vector2<double> disp = other_obj.pos - this_obj.pos;
				
				const double dist_sq = disp.x * disp.x + disp.y * disp.y;
				const double dist = sqrt(dist_sq);
				const double inv_dist_sq = 1 / dist_sq;
				const double inv_dist = 1 / dist;
				const sf::Vector2<double> disp_norm = disp / dist;
				const sf::Vector2<double> theta_dir = {disp_norm.y, -disp_norm.x};

				const sf::Vector2<double> other_dir = { sin(other_obj.angle), cos(other_obj.angle) };

				const double this_dot_r = disp_norm.x * this_dir.x + disp_norm.y * this_dir.y;
				const double other_dot_r = disp_norm.x * other_dir.x + disp_norm.y * other_dir.y;
				const double this_dot_other = this_dir.x * other_dir.x + this_dir.y * other_dir.y;

				const double this_cross_other = this_dir.x * other_dir.y - this_dir.y * other_dir.x;
				const double this_cross_r = this_dir.x * disp_norm.y - this_dir.y * disp_norm.x;
				const double other_cross_r = other_dir.x * disp_norm.y - other_dir.y * disp_norm.x;

				const double r_component = this_dot_other - 3 * this_dot_r * other_dot_r;
				const double theta_component = this_dot_r * other_cross_r + other_dot_r * this_cross_r; // these aren't actually cross products, but since they're already computed it's simpler to do it this way.

				//std::cout << r_component << ", " << theta_component << std::endl;

				const double force_factor = dipole_k * inv_dist_sq * inv_dist_sq;
				const double torque_factor = -dipole_k * .5 * inv_dist_sq * inv_dist;

				//std::cout << force_factor << ", " << factor << std::endl;

				const sf::Vector2<double> force = force_factor * (r_component * disp_norm + theta_component * theta_dir);
				const double torque_on_other = torque_factor * (3 * this_dot_r * other_cross_r + this_cross_other);
				const double torque_on_this = torque_factor * (3 * other_dot_r * this_cross_r - this_cross_other);

				//other_obj.addAcc(force); // assuming the dipole moment is 1 for all objects for now
				//this_obj.addAcc(-force);

				other_obj.addAngAcc(torque_on_other);
				this_obj.addAngAcc(torque_on_this);
			}
		}
	}

	void updateObjects(double dt) {
		for (PhysicsObject& obj : objects) {
			obj.update(dt);
		}
	}
};