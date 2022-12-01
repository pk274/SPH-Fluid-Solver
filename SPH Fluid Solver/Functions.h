// Paul Kull, 2022
#pragma once

#include<SFML/Graphics.hpp>


class Functions {
 public:
	static sf::Vector2f calculate_distance(sf::Vector2f pos1, sf::Vector2f pos2);
	static double calculate_distance_norm(sf::Vector2f distance);
	static double kernel(double distance);
	static sf::Vector2f kernel_derivation(sf::Vector2f distance, float distanceNorm);
	static double scalar_product2D(sf::Vector2f v1, sf::Vector2f v2);
	static float round(double number, int places = 3);
};