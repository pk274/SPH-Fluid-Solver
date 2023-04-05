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
	  sf::Vector2f _velocity;
	  sf::Vector2f _pressureAcc;

	  sf::Vector2f _v_adv;
	  sf::Vector2f _x_adv;
	  sf::Vector2f c_f;

	  float _density;
	  float _pressure;
	  int _colorFactor;

	  float _a_ii;
	  float _s_i;

	  std::vector<Particle*> _neighbors;

	  // Rendering Information
	  sf::Color _stasisColor;

};