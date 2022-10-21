// Paul Kull, 2022
#pragma once

#include <SFML/Graphics.hpp>

enum ParticleType
{
	solid = 0,
	fluid = 1
};

class Particle {
  public:
	  // Particle Type
	  ParticleType _type;

	  // Physical Information
	  sf::Vector2f _position;

	  // Rendering Information
	  sf::Color _stasisColor;

};