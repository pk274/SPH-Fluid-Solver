// Paul Kull, 2022

#define _USE_MATH_DEFINES

#include <iostream>
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
std::vector<int> TestManager::test_kernel(std::vector<Particle>* particles, float h, int watchedId) {
	std::vector<int> positivelyTestedParticles = std::vector<int>();
	int numNeighbors;
	sf::Vector2f kernelDerivativeSum;
	float areaSum[2][2];
	sf::Vector2f distance_ij;
	double distance_ijNorm;
	sf::Vector2f distance_ji;
	double distance_jiNorm;
	sf::Vector2f kernelDeriv;
	bool symmetric = true;
	bool sumZero;
	bool areaIsRight;
	float area = 7.32600732601;
	for (int i = 0; i < particles->size(); i++) {
		// if (particles->at(i)._type == solid) {continue; }
		kernelDerivativeSum.x = 0;
		kernelDerivativeSum.y = 0;
		areaSum[0][0] = 0;
		areaSum[0][1] = 0;
		areaSum[1][0] = 0;
		areaSum[1][1] = 0;
		numNeighbors = particles->at(i)._neighbors.size();
		for (int j = 0; j < numNeighbors; j++) {
			distance_ij = Functions::calculate_distance(particles->at(i)._position, particles->at(i)._neighbors.at(j)->_position);
			distance_ijNorm = Functions::calculate_distance_norm(distance_ij);
			distance_ji = Functions::calculate_distance(particles->at(i)._neighbors.at(j)->_position, particles->at(i)._position);
			distance_jiNorm = Functions::calculate_distance_norm(distance_ji);
			kernelDeriv = Functions::kernel_derivation(distance_ij, distance_ijNorm, h);

			if (!(kernelDeriv == -Functions::kernel_derivation(distance_ji, distance_jiNorm, h))) {
				symmetric = false;
			}
			kernelDerivativeSum += kernelDeriv;
			areaSum[0][0] += distance_ij.x * kernelDeriv.x * 1 / numNeighbors;
			areaSum[0][1] += distance_ij.x * kernelDeriv.y * 1 / numNeighbors;
			areaSum[1][0] += distance_ij.y * kernelDeriv.x * 1 / numNeighbors;
			areaSum[1][1] += distance_ij.y * kernelDeriv.y * 1 / numNeighbors;
		}
		sumZero = (abs(kernelDerivativeSum.x) < THRESHHOLD) && (abs(kernelDerivativeSum.y) < THRESHHOLD);
		areaIsRight = true;
		// This stupid fuk grows with growing kernel area
		areaIsRight = ((abs(areaSum[0][0] + 1 / area) < THRESHHOLD) && (abs(areaSum[1][1] + 1 / area) < THRESHHOLD)
			&& (abs(areaSum[0][1] < THRESHHOLD)) && (abs(areaSum[1][0] < THRESHHOLD)));
		if (symmetric && sumZero && areaIsRight) {
			positivelyTestedParticles.push_back(particles->at(i)._id);
		}
		// if (i == watchedId) {
		// 	std::cout << areaSum[0][0] << "	" << areaSum[0][1] << std::endl;
		// 	std::cout << areaSum[1][0] << "	" << areaSum[1][1] << "\n" << std::endl;
		// }
	}
	return positivelyTestedParticles;
}

// _________________________________________________________________________________________
void TestManager::test_kernel_integral() {
	int areaSize = 6;
	float stepSize = 2;
	float H = 2;
	double kernelIntegral = 0;
	double kernelSum = 0;
	for (float i = -areaSize; i <= areaSize; i += stepSize) {
		for (float j = -areaSize; j <= areaSize; j += stepSize) {
			if (i == 0 && j == 0) { continue; }
			kernelIntegral += stepSize * stepSize * Functions::kernel(
				Functions::calculate_distance_norm(sf::Vector2f(i, j)), H);
			kernelSum += Functions::kernel(
				Functions::calculate_distance_norm(sf::Vector2f(i, j)), H);
			std::cout << 4 * Functions::kernel(
				Functions::calculate_distance_norm(sf::Vector2f(i, j)), H) << std::endl;
		}
	}
	std::cout << "\n" << kernelIntegral << std::endl;
	std::cout << 4 * kernelSum << std::endl;
}
