#include <chrono>
#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include "renderer.h"
#include "solver.h"

void initRandom(Solver& solver, int lowx, int lowy, int highx, int highy, int num_objects, double radius) {
	unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_int_distribution<int> distribution_x(lowx, highx);
	std::uniform_int_distribution<int> distribution_y(lowy, highy);

	

	for (int i = 0; i < num_objects; i++) {
		double x = distribution_x(generator);
		double y = distribution_y(generator);
		//std::cout << x << ", " << y << std::endl;
		//solver.addObject(PhysicsObject({x, y}, 5));
		solver.addObject({ x, y }, 0, radius);
	}
}

void initSquare(Solver& solver, int lowx, int lowy, int highx, int highy, int num_objects, double radius, int distx, double disty, double rand_amt, bool angle) {
	unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution_x(-rand_amt, rand_amt);
	std::uniform_real_distribution<double> distribution_y(-rand_amt, rand_amt);
	std::uniform_real_distribution<double> rand_angle(0, 6.28318530718);

	for (int i = 0; i < num_objects; i++) {
		double dx = highx - lowx;
		double dy = highy - lowy;
		int j = i;
		double x = distribution_x(generator) + lowx + (distx * j) % (int)dx;
		double y = distribution_y(generator) + lowy + disty * floor(distx * j / dy);
		double a = 0;
		if (angle) a = rand_angle(generator);
		//double x = lowx + i * dx / num_objects;
		//double y = lowy + floor(i / dy);
		//std::cout << x << ", " << y << std::endl;
		//solver.addObject(PhysicsObject({x, y}, 5));
		solver.addObject({ x, y }, a, radius);
	}
}

void initSquareNew(Solver& solver, double start_x, double start_y, int row_size, int num_objects, double radius, double distx, double disty, double rand_amt, bool angle, bool charged) {
	unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution_x(-rand_amt, rand_amt);
	std::uniform_real_distribution<double> distribution_y(-rand_amt, rand_amt);
	std::uniform_real_distribution<double> rand_angle(0, 6.28318530718);

	for (int i = 0; i < num_objects; i++) {
		int j = i;
		double x = distribution_x(generator) + start_x + distx * (j % row_size);
		double y = distribution_y(generator) + start_y + disty * floor(j / row_size);
		double a = 0;
		if (angle) a = rand_angle(generator);
		//double x = lowx + i * dx / num_objects;
		//double y = lowy + floor(i / dy);
		//std::cout << x << ", " << y << std::endl;
		//solver.addObject(PhysicsObject({x, y}, 5));
		double charge = 0;
		if (charged) charge = 1;//floor((i % 3) / 2) * 2 - 1;// (i % 2) * 2 - 1;
		solver.addObject({ x, y }, a, charge, radius);
	}
}

bool contains(sf::RectangleShape& shape, sf::Vector2i& p) {
	return (p.x >= shape.getPosition().x && p.x <= shape.getPosition().x + shape.getSize().x && p.y >= shape.getPosition().y && p.y <= shape.getPosition().y + shape.getSize().y);
}

int startScreen(sf::RenderWindow& window, Solver& solver) {
	sf::RectangleShape start_button({ 200, 50 });
	start_button.setFillColor(sf::Color::White);
	start_button.setPosition({window.getSize().x - start_button.getSize().x - .5f * start_button.getSize().y, window.getSize().y - 1.5f * start_button.getSize().y });

	sf::Font font;
	font.loadFromFile("fonts/ARIAL.ttf");
	sf::Text start_button_text("Start", font, 30);
	start_button_text.setFillColor(sf::Color::Black);
	start_button_text.setPosition({start_button.getPosition().x + 50, start_button.getPosition().y + 5});

	bool do_dipole = false;
	solver.setDoDipole(do_dipole);
	sf::RectangleShape dipole_button({400, 50});
	dipole_button.setPosition({50, 50});
	sf::Text dipole_button_text("Turn on dipoles? green = true", font, 20);
	dipole_button_text.setFillColor(sf::Color::Black);
	dipole_button_text.setPosition({ dipole_button.getPosition().x, dipole_button.getPosition().y + 5 });

	/*
	double lj_r0 = 20;
	solver.setLJR0(lj_r0);
	sf::RectangleShape lj_r0_button({ 400, 50 });
	lj_r0_button.setPosition({ 50, 150 });
	lj_r0_button.setFillColor(sf::Color::White);
	sf::Text lj_r0_button_text("set LJ equilibrium radius (default = 20):", font, 20);
	lj_r0_button_text.setFillColor(sf::Color::Black);
	lj_r0_button_text.setPosition({ lj_r0_button.getPosition().x, lj_r0_button.getPosition().y + 5 });
	*/

	int num_particles = 200;
	sf::RectangleShape num_particles_button({ 400, 50 });
	num_particles_button.setPosition({ 50, 150 });
	num_particles_button.setFillColor(sf::Color::White);
	sf::Text num_particles_button_text("num particles? white = 200, green = 500:", font, 20);
	num_particles_button_text.setFillColor(sf::Color::Black);
	num_particles_button_text.setPosition({ num_particles_button.getPosition().x, num_particles_button.getPosition().y + 5 });

	//std::cout << window.getSize().y - 1.5f * start_button.getSize().y << std::endl;
	//std::cout << start_button_text.getFont() << std::endl;
	while (window.isOpen()) {
		sf::Event event;
		while (window.pollEvent(event)) {
			if (event.type == sf::Event::Closed || sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
				window.close();
			}
			if (event.type == sf::Event::MouseButtonPressed) {
				sf::Vector2i mouse_pos = sf::Mouse::getPosition(window);
				if (contains(start_button, mouse_pos)) {
					return num_particles; // this is the most jank thing ever "i will fix later" - me
				}

				if (contains(dipole_button, mouse_pos)) {
					do_dipole = (do_dipole + 1) % 2;
					//std::cout << do_dipole << std::endl;
					solver.setDoDipole(do_dipole);
				} 
				if (contains(num_particles_button, mouse_pos)) {
					if (num_particles == 200) {
						num_particles = 500;
					} else {
						num_particles = 200;
					}
				}
			}
		}
		window.clear(sf::Color::Color(100, 100, 100));
		window.draw(start_button);
		window.draw(start_button_text);

		if (do_dipole) {
			dipole_button.setFillColor(sf::Color::Green);
		} else {
			dipole_button.setFillColor(sf::Color::White);
		}
		window.draw(dipole_button);
		window.draw(dipole_button_text);

		if (num_particles == 500) {
			num_particles_button.setFillColor(sf::Color::Green);
		} else {
			num_particles_button.setFillColor(sf::Color::White);
		}
		window.draw(num_particles_button);
		window.draw(num_particles_button_text);

		//window.draw(lj_r0_button);
		//window.draw(lj_r0_button_text);

		window.display();


	}
}

int main() {
	const double pi = 3.14159265358979323846; // this level of precision is pretty dumb
	const int width = 800;
	const int height = 800;
	const double r = 10;
	bool do_capture = false;
	bool do_dipole = false;

	//const int num_cols = width / cell_size;
	//const int num_rows = height / cell_size;

	//int grid[num_rows][num_cols];

	//for (int i = 0; i < num_rows; i++) {
	//	for (int j = 0; j < num_cols; j++) {
	//		if ((i * j + i + j) % 5 == 0) {
	//			grid[i][j] = 1;
	//		} else {
	//			grid[i][j] = 0;
	//		}
	//	}
	//}
	//std::cout << std::chrono::system_clock::now().time_since_epoch().count() << std::endl;
	//std::cout << std::chrono::system_clock::now().time_since_epoch().count() << std::endl;
	sf::ContextSettings settings;
	settings.antialiasingLevel = 4;

	sf::RenderWindow window(sf::VideoMode(width, height), "bruh", sf::Style::Default, settings);

	const int framerate = 60;
	//window.setFramerateLimit(framerate);



	Solver solver(window.getSize(), do_dipole);
	Renderer renderer{ window };
	const int subs_steps = 10;
	solver.setSubStepsCount(subs_steps);
	solver.setSimUpdateRate(framerate);

	int num_particles = startScreen(window, solver);
	do_dipole = solver.getDoDipole();
	// this do_dipole stuff is the most serious care of malarkey i've ever seen

	//PhysicsObject test1 = PhysicsObject({ 220, height / 2 }, 5);
	//test1.addVel({0, -10}, 1. / (subs_steps * framerate));

	//PhysicsObject test1 = PhysicsObject({ 205, height / 2 }, 2);
	//PhysicsObject test2 = PhysicsObject({ 195, height / 2 }, 2);

	//solver.addObject({ width * .5 + 25, height / 2 }, .5 * pi, 10);
	//solver.addObject({ width * .5 - 25, height / 2 }, 0, 10);
	//double start_energy = solver.calcTotalEnergyLJ_debug();
	//std::cout << start_energy << std::endl;

	//initRandom(solver, width * .2, height * .2, width * .8, height * .8, 50, r);
	//initSquare(solver, width * .15, height * .15, width * .85, height * .85, 500, r, 25, 25, 5);
	//initSquare(solver, width * .2, height * .2, width * .8, height * .8, 300, r, 25, 21, 0, true);
	//initSquareNew(solver, width * .1, height * .2, 25, 500, r, 11, 11, 1, true, true);
	initSquareNew(solver, width * .05, height * .2, 38, num_particles, r, 20, 20, 1, do_dipole, false);


	int counter = 0;
	int frame_count = 0;
	bool outin = true;

	//std::cout << solver.getStepDT() << std::endl;

	auto tp1 = std::chrono::high_resolution_clock::now();

	//std::cout << )((std::chrono::system_clock::now() - tp1)) << std::endl;

	while (window.isOpen()) {
		tp1 = std::chrono::high_resolution_clock::now();
		sf::Event event;
		while (window.pollEvent(event)) {
			if (event.type == sf::Event::Closed || sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
				window.close();
			}
			if (event.type == sf::Event::MouseButtonPressed) {
				unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
				std::default_random_engine generator(seed);
				std::uniform_real_distribution<double> rand_angle(0, 6.28318530718);
				sf::Vector2<double> m_pos = static_cast<sf::Vector2<double>>(sf::Mouse::getPosition(window));
				solver.addObject({ m_pos.x, m_pos.y }, rand_angle(generator), r);
			}
		}

		solver.update();
		window.clear(sf::Color::Color(100, 100, 100));
		renderer.render(solver, do_dipole);
		window.display();

		//const PhysicsObject& test1_ref = solver.getObjects()[0];
		//const PhysicsObject& test2_ref = solver.getObjects()[1];
		
		//std::cout << test2_ref.pos.x << ", " << test2_ref.pos.y << std::endl;
		//std::cout << test1_ref.angle << std::endl;
		//double temp_e = solver.calcTotalEnergyLJ_debug();
		//std::cout << "AHHHHHHHHH " << temp_e << ", " << counter << std::endl;
		//if (temp_e != start_energy) {
			//std::cout << "AHHHHHHHHH " << temp_e << ", " << counter << std::endl;
		//}
		//sf::Vector2<double> v = test1_ref.getVel(1. / (subs_steps * framerate));
		//std::cout << test1_ref.pos.x << ", " << v.x << ", " << temp_e << std::endl;
		//std::cout << test1_ref.pos.x << ", " << test1_ref.pos.y  << ", " << temp_e << std::endl;
		//std::cout << test1_ref.pos.y << std::endl;

		/*
		if (test1_ref.pos.x > width * .5 + 25 && outin) {
			counter++;
			std::cout << "AHHHHHHHHH " << test1_ref.pos.x << ", " << counter << std::endl;
			//outin = false;
		}
		if (!(test1_ref.pos.x > width * .5 + 25)) {
			outin = true;
		}
		*/
		//std::cout << "frame count : " << frame_count << std::endl;
		
		if (do_capture) {
			sf::Vector2u window_size = window.getSize();
			sf::Texture texture;
			texture.create(window_size.x, window_size.y);
			texture.update(window);
			sf::Image screenshot = texture.copyToImage();
			//screenshot.saveToFile("pics/test" + std::to_string(frame_count) + ".png");
		}
		
		frame_count++;
		double millis = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - tp1).count() / 1000000.;
		//std::cout << millis << "ms per frame, " << 1000 / millis << " fps, " << frame_count << " frames" << std::endl;
		std::cout << millis << "ms per frame, " << 1000/millis << " fps, " << frame_count << " frames, " << solver.getTotalVelSq(1. / (subs_steps * framerate)) / solver.getObjectsCount() <<  "avg total vel squared" << std::endl; // frame time in millis
	}
	return 0;
}