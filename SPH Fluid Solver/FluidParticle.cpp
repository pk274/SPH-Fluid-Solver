// Paul Kull, 2022

#include "./FluidParticle.h"


// _________________________________________________________________
FluidParticle::FluidParticle(int id, sf::Vector2f pos) {
	// Initialize Type
	_type = fluid;
	_id = id;

	// Initialize Position
	_position.x = pos.x;
	_position.y = pos.y;


	// Init Rendering Information
	_stasisColor = sf::Color::Blue;


}