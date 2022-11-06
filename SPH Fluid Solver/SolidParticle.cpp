// Paul Kull, 2022

#include "./SolidParticle.h"


const float SolidParticle::_size = 2.f;
const float SolidParticle::_mass = 8.f;


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


	// Init Rendering Information
	_stasisColor = sf::Color::White;


}