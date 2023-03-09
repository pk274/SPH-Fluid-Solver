// Paul Kull, 2022

#include "./FluidParticle.h"


const float FluidParticle::_size = 2.f;	// cm Durchmesser
const float FluidParticle::_mass = 4.f;		// g
const float FluidParticle::_restDensity = 1.f;
const float FluidParticle::_materialParameter = 1.f;

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

	_v_adv = sf::Vector2f();


	// Init Rendering Information
	_stasisColor = sf::Color::Blue;

}