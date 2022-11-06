// Paul Kull, 2022

#define _USE_MATH_DEFINES

#include "./TestManager.h"

constexpr double THRESHHOLD = 0.00001;


// ______________________________________________________________________________
std::vector<int> TestManager::test_correct_neighbor_amount(std::vector<Particle>* particles, int supposedNumNeighbors) {
	std::vector<int> positivelyTestedParticles = std::vector<int>();
	for (int i = 0; i < particles->size(); i++) {
		if (particles->at(i)._neighbors.size() == supposedNumNeighbors) {
			positivelyTestedParticles.push_back(particles->at(i)._id);
		}
	}
	return positivelyTestedParticles;
}

// ______________________________________________________________________________
std::vector<int> TestManager::test_kernel(std::vector<Particle>* particles, float h) {
	std::vector<int> positivelyTestedParticles = std::vector<int>();
	sf::Vector2f kernelDerivativeSum;
	sf::Vector2f areaSum;
	sf::Vector2f distance_ij;
	double distance_ijNorm;
	sf::Vector2f distance_ji;
	double distance_jiNorm;
	sf::Vector2f kernelDeriv;
	bool symmetric = true;
	bool sumZero;
	bool areaIsRight;
	float area = 4;
	for (int i = 0; i < particles->size(); i++) {
		// if (particles->at(i)._type == solid) {continue; }
		kernelDerivativeSum.x = 0;
		kernelDerivativeSum.y = 0;
		areaSum.x = 0;
		areaSum.y = 0;
		for (int j = 0; j < particles->at(i)._neighbors.size(); j++) {
			distance_ij = Functions::calculate_distance(particles->at(i)._position, particles->at(i)._neighbors.at(j)->_position);
			distance_ijNorm = Functions::calculate_distance_norm(distance_ij);
			distance_ji = Functions::calculate_distance(particles->at(i)._neighbors.at(j)->_position, particles->at(i)._position);
			distance_jiNorm = Functions::calculate_distance_norm(distance_ji);
			kernelDeriv = Functions::kernel_derivation(distance_ij, distance_ijNorm, h);

			if (!(kernelDeriv == - Functions::kernel_derivation(distance_ji, distance_jiNorm, h))) {
				// Usually not a problem
				symmetric = false;
			}
			kernelDerivativeSum += kernelDeriv;
			areaSum.x += distance_ij.x * kernelDeriv.x;
			areaSum.y += distance_ij.y * kernelDeriv.y;
		}
		sumZero = (abs(kernelDerivativeSum.x) < THRESHHOLD) && (abs(kernelDerivativeSum.y) < THRESHHOLD);
		areaIsRight = true;
		areaIsRight = ((abs(areaSum.x + 1 / area) < THRESHHOLD) && (abs(areaSum.y + 1 / area) < THRESHHOLD));
		if (symmetric && sumZero && areaIsRight) {
			positivelyTestedParticles.push_back(particles->at(i)._id);
		}
	}
	return positivelyTestedParticles;
}
