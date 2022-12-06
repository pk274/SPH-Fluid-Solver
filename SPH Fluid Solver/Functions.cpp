// Paul Kull, 2022

#define _USE_MATH_DEFINES

#include "./Functions.h"
#include "./Parameters.h"

constexpr float kernelCorrection = 0.99914073896;
float a = 5 / (14 * M_PI * Parameters::H * Parameters::H);


// __________________________________________________________________________________
sf::Vector2f Functions::calculate_distance(sf::Vector2f pos1, sf::Vector2f pos2) {
	return sf::Vector2f(pos1.x - pos2.x, pos1.y - pos2.y);
}


// __________________________________________________________________________________
double Functions::calculate_distance_norm(sf::Vector2f distance) {
	return std::sqrt(std::pow(distance.x, 2) + std::pow(distance.y, 2));
}



// _________________________________________________________________________________
double Functions::kernel(double distance) {
	double q = distance / Parameters::H;
	double t1 = std::max(1 - q, 0.);
	double t2 = std::max(2 - q, 0.);
	return a * kernelCorrection * (t2 * t2 * t2 - 4 * t1 * t1 * t1);
}

// _________________________________________________________________________________
sf::Vector2f Functions::kernel_derivation(sf::Vector2f distance, float distanceNorm) {
	if (distanceNorm == 0) { return sf::Vector2f(0, 0); }
	float q = distanceNorm / Parameters::H;
	float t1 = std::max(1 - q, 0.f);
	float t2 = std::max(2 - q, 0.f);
	return a * distance / (distanceNorm * Parameters::H) * (-3 * t2 * t2 + 12 * t1 * t1);
}

// _________________________________________________________________________________
double Functions::scalar_product2D(sf::Vector2f v1, sf::Vector2f v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

// _________________________________________________________________________________
float Functions::round(double number, int places) {
	int factor = std::pow(10, places);
	return std::round(factor * number) / factor;
}
