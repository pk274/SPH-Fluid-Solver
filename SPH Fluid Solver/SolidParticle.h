// Paul Kull, 2022
#pragma once

#include <SFML/Graphics.hpp>

#include "./Particle.h"

class SolidParticle : public Particle {
  public:
	static const int _size = 50;

	SolidParticle(sf::Vector2f pos);
};