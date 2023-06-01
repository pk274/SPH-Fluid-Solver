// Paul Kull, 2023

#pragma once

#include <SFML/Graphics.hpp>
#include "./Particle.h"

class MovingObject {
public:
	std::vector<sf::Vector2f> _conditions;
	std::vector<bool> _conditionBigger;
	std::vector<sf::Vector2f> _directions;
	std::vector<Particle*> _particles;
	int _state;
};