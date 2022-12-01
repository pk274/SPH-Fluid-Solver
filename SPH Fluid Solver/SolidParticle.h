// Paul Kull, 2022
#pragma once

#include <SFML/Graphics.hpp>

#include "./Particle.h"

class SolidParticle : public Particle {
  public:
	static const float _size;
	static const float _mass;

	// Rendering Information
	const static sf::Color _stasisColor;

	SolidParticle(int id, sf::Vector2f pos);
};

