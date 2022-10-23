// Paul Kull, 2022

#include "./SolidParticle.h"


// _________________________________________________________________
SolidParticle::SolidParticle(int id, sf::Vector2f pos) {
	// Initialize Type
	_type = solid;
	_id = id;

	// Initialize Position
	_position.x = pos.x;
	_position.y = pos.y;

	// Init Rendering Information
	_stasisColor = sf::Color::White;


}