// Paul Kull, 2022
#pragma once

#include<SFML/Graphics.hpp>
#include "./Particle.h"


class Functions {
 public:
	static sf::Vector2f calculate_distance(sf::Vector2f pos1, sf::Vector2f pos2);
	static float calculate_distance_norm(sf::Vector2f distance);
	static float kernel(float distance);
	static sf::Vector2f kernel_derivation(sf::Vector2f distance, float distanceNorm);
	static float scalar_product2D(sf::Vector2f v1, sf::Vector2f v2);
	static float round(float number, int places = 3);
	static std::vector<Particle*> n_square_neighborhood_search(std::vector<Particle>* particles, int index, float radius);
	static std::tuple<int, int, int> color_code_pressure(float pressure);
};