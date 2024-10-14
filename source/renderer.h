#pragma once
#include <SFML/Graphics.hpp>
#include <iostream>
#include "solver.h"

class Renderer {
	const double pi = 3.14159265358979323846;
	const double to_degrees = 180 / pi;
public:
	sf::Vector2f window_size;
	explicit Renderer(sf::RenderTarget& target) : render_target{target} {
		window_size = static_cast<sf::Vector2f>(render_target.getSize());
	}

	void render(const Solver& solver, bool do_dipole) const {

		sf::CircleShape circle(1.f);
		circle.setPointCount(16);
		circle.setOrigin(1.f, 1.f);
		
		
		circle.setPosition(window_size * .5f);
		circle.setScale(2, 2);
		circle.setFillColor(sf::Color::Magenta);
		render_target.draw(circle);

		sf::ConvexShape cs;
		cs.setPointCount(3);
		cs.setPoint(0, { 0.f, -1.f });
		cs.setPoint(1, { -.55f, -.25f });
		cs.setPoint(2, { .55f, -.25f });

		sf::RectangleShape line({cs.getPoint(2).x * .6f, 1.22});
		line.setOrigin({line.getSize().x * .5f, -cs.getPoint(2).y});

		const std::vector<PhysicsObject>& objects = solver.getObjects();
		for (const PhysicsObject& obj : objects) {

			circle.setPosition(static_cast<sf::Vector2f>(obj.pos));
			circle.setScale(obj.radius, obj.radius);
			if (obj.charge == 0) {
				if (do_dipole) {
					circle.setFillColor(obj.getColorDipole());
				} else {
					//circle.setFillColor(sf::Color::White);
					const double obj_vel_mag = obj.getVelMagSq(solver.getStepDT());
					const double clamped = std::max(255 - obj_vel_mag, 0.);
					//circle.setFillColor(hsv(obj_vel_mag * 1 + 180, 1, 1));
					circle.setFillColor(sf::Color::Color(255, clamped, clamped));
				}
			} else {
				circle.setFillColor(obj.getColorCharge());
			}
			render_target.draw(circle);
			
			// this is for the dipole
			if (do_dipole) {
				if (obj.charge == 0) {
					cs.setPosition(static_cast<sf::Vector2f>(obj.pos));
					cs.setScale(obj.radius, obj.radius);
					cs.setRotation(obj.angle * to_degrees);
					cs.setFillColor(sf::Color::Black);
					render_target.draw(cs);

					line.setPosition(static_cast<sf::Vector2f>(obj.pos));
					line.setScale(obj.radius, obj.radius);
					line.setRotation(obj.angle * to_degrees);
					line.setFillColor(sf::Color::Black);
					render_target.draw(line);
				}
			}
			
		}
	}


private:
	sf::RenderTarget& render_target;
};
