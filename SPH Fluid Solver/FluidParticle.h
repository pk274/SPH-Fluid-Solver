// Paul Kull, 2022
#pragma once

#include <SFML/Graphics.hpp>

#include "./Particle.h"


class FluidParticle : public Particle {
  public:
	static const int _size = 1;

	FluidParticle(int id, sf::Vector2f pos);
};
