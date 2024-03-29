// Paul Kull, 2022

#define _USE_MATH_DEFINES

#include "./Functions.h"
#include "./Parameters.h"

constexpr float kernelCorrection = 0.99914073896;
// constexpr float kernelCorrection = 1;
float a = 5 / (14 * M_PI * Parameters::H * Parameters::H);


// __________________________________________________________________________________
sf::Vector2f Functions::calculate_distance(sf::Vector2f pos1, sf::Vector2f pos2) {
	return sf::Vector2f(pos1.x - pos2.x, pos1.y - pos2.y);
}


// __________________________________________________________________________________
float Functions::calculate_distance_norm(sf::Vector2f distance) {
	return std::sqrt(distance.x * distance.x + distance.y * distance.y);
}

// This implementation is largely taken from the slides of the course on fluid simulation
// by the computer graphics department of the university of Freiburg.
// _________________________________________________________________________________
float Functions::kernel(float distance) {
	float q = distance / Parameters::H;
	float t1 = std::max(1 - q, 0.f);
	float t2 = std::max(2 - q, 0.f);
	return a * kernelCorrection * (t2 * t2 * t2 - 4 * t1 * t1 * t1);
}

// This implementation is largely taken from the slides of the course on fluid simulation
// by the computer graphics department of the university of Freiburg.
// _________________________________________________________________________________
sf::Vector2f Functions::kernel_derivation(sf::Vector2f distance, float distanceNorm) {
	if (distanceNorm == 0) { return sf::Vector2f(0, 0); }
	float q = distanceNorm / Parameters::H;
	float t1 = std::max(1 - q, 0.f);
	float t2 = std::max(2 - q, 0.f);
	return a * distance / (distanceNorm * Parameters::H) * (-3 * t2 * t2 + 12 * t1 * t1);
}

// _________________________________________________________________________________
float Functions::scalar_product2D(sf::Vector2f v1, sf::Vector2f v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

// _________________________________________________________________________________
float Functions::round(float number, int places) {
	int factor = std::pow(10, places);
	return std::round(factor * number) / factor;
}


// _________________________________________________________________________________
std::vector<Particle*> Functions::n_square_neighborhood_search(std::vector<Particle>* particles, int index, float radius) {
	std::vector<Particle*> neighbors = std::vector<Particle*>();
	for (int i = 0; i < particles->size(); i++) {
		if (calculate_distance_norm(calculate_distance
		(particles->at(i)._position, particles->at(index)._position)) < radius) {
			neighbors.push_back(&particles->at(i));
		}
	}
	return neighbors;
}

// ___________________________________________________________________________________________
std::tuple<int, int, int> Functions::color_code_pressure(float pressure) {
	if (pressure <= 1) {
		return std::make_tuple<int, int, int>(48, 33, 217);
	}
	else if (pressure <= 300 * Parameters::PRESSURE_CODE_ROUGHNESS) {
		return std::make_tuple<int, int, int>(33, 90, 217);
	}
	else if (pressure <= 800 * Parameters::PRESSURE_CODE_ROUGHNESS) {
		return std::make_tuple<int, int, int>(33, 130, 217);
	}
	else if (pressure <= 2000 * Parameters::PRESSURE_CODE_ROUGHNESS) {
		return std::make_tuple<int, int, int>(33, 162, 217);
	}
	else if (pressure <= 5000 * Parameters::PRESSURE_CODE_ROUGHNESS) {
		return std::make_tuple<int, int, int>(33, 217, 186);
	}
	else if (pressure <= 10000 * Parameters::PRESSURE_CODE_ROUGHNESS) {
		return std::make_tuple<int, int, int>(33, 217, 125);
	}
	else if (pressure <= 25000 * Parameters::PRESSURE_CODE_ROUGHNESS) {
		return std::make_tuple<int, int, int>(33, 217, 42);
	}
	else if (pressure <= 50000 * Parameters::PRESSURE_CODE_ROUGHNESS) {
		return std::make_tuple<int, int, int>(131, 217, 33);
	}
	else if (pressure <= 70000 * Parameters::PRESSURE_CODE_ROUGHNESS) {
		return std::make_tuple<int, int, int>(199, 217, 33);
	}
	else if (pressure <= 90000 * Parameters::PRESSURE_CODE_ROUGHNESS) {
		return std::make_tuple<int, int, int>(217, 174, 33);
	}
	else if (pressure <= 110000 * Parameters::PRESSURE_CODE_ROUGHNESS) {
		return std::make_tuple<int, int, int>(217, 125, 33);
	}
	else if (pressure <= 140000 * Parameters::PRESSURE_CODE_ROUGHNESS) {
		return std::make_tuple<int, int, int>(217, 79, 33);
	}
	else if (pressure <= 170000 * Parameters::PRESSURE_CODE_ROUGHNESS) {
		return std::make_tuple<int, int, int>(217, 48, 33);
	}
	else {
		return std::make_tuple<int, int, int>(181, 18, 4);
	}
}