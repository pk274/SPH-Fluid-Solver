// Paul Kull, 2022

#include "./SolidParticle.h"

// _________________________________________________________________
SolidParticle::SolidParticle(sf::Vector2f pos) {
	// Initialize Type
	_type = solid;

	// Initialize Position
	_position.x = pos.x;
	_position.y = pos.y;

	// Init Rendering Information
	_stasisColor = sf::Color::White;


}