// Paul Kull, 2022

#include "./FluidParticle.h"


const float FluidParticle::_size = 1.;
const float FluidParticle::_mass = 1.;
const float FluidParticle::_restDensity = 1.;

// _________________________________________________________________
FluidParticle::FluidParticle(int id, sf::Vector2f pos) {
	// Initialize Type
	_type = fluid;
	_id = id;

	// Initialize Position
	_position.x = pos.x;
	_position.y = pos.y;
	_velocity.x = 0;
	_velocity.x = 0;


	// Init Rendering Information
	_stasisColor = sf::Color::Blue;


}