// Paul Kull, 2022

#include "./SolidParticle.h"


const float SolidParticle::_size = 2.f;
const float SolidParticle::_mass = 4.f;


// _________________________________________________________________
SolidParticle::SolidParticle(int id, sf::Vector2f pos) {
	// Initialize Type
	_type = solid;
	_id = id;

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