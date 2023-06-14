// Paul Kull 2023

#include "./MovingSolidParticle.h"


const float MovingSolidParticle::_size = 2.f;
const float MovingSolidParticle::_mass = 4.f;


// _________________________________________________________________
MovingSolidParticle::MovingSolidParticle(int id, sf::Vector2f pos, sf::Vector2f vel, std::vector<sf::Vector2f> velocities,
	std::vector<sf::Vector2f> triggers, std::vector<int> friends) {
	// Initialize Type
	_type = movingSolid;
	_id = id;
	_statusPointer = 0;

	// Initialize Position
	_position.x = pos.x;
	_position.y = pos.y;

	_velocity = sf::Vector2f(0, 0);

	_density = 1;
	_pressure = 0;

	_v_adv = sf::Vector2f(0, 0);
	_x_adv = _position;
	_pressureAcc = sf::Vector2f(0, 0);
	c_f = sf::Vector2f(0, 0);

	_a_ii = 0;
	_s_i = 0;

	_colorFactor = 0;

	// Init render info
	_stasisColor = sf::Color::White;

}