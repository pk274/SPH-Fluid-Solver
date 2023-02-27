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

	_velocity.x = 0;
	_velocity.x = 0;

	_density = 1;
	_pressure = 0;

	_colorFactor = 0;

	// Init render info
	_stasisColor = sf::Color::White;

}