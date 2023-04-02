// Paul Kull, 2022
#pragma once

#include <SFML/Graphics.hpp>

#include "./Particle.h"


class FluidParticle : public Particle {
  public:
	static const float _size;
	static const float _mass;
	static const float _restDensity;
	static const float _materialParameter;

	FluidParticle(int id, sf::Vector2f pos, sf::Vector2f vel = sf::Vector2f(0, 0));
};
