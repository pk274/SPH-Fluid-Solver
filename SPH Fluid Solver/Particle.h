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
	  int _id;

	  // Physical Information
	  sf::Vector2f _position;

	  int _distance;

	  // Rendering Information
	  sf::Color _stasisColor;

	  sf::Time _lastUpdated;

};