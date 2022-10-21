// Paul Kull, 2022

#include "./FluidParticle.h"

// _________________________________________________________________
FluidParticle::FluidParticle(sf::Vector2f pos) {
	// Initialize Type
	_type = fluid;

	// Initialize Position
	_position.x = pos.x;
	_position.y = pos.y;

	// Init Rendering Information
	_stasisColor = sf::Color::Blue;


}