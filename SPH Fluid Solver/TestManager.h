// Paul Kull, 2022

#pragma once

#include <vector>
#include "./TestManager.h"
#include "./Functions.h"
#include "./Particle.h"

class TestManager {
  public:
	static std::vector<int> test_correct_neighbor_amount(std::vector<Particle>* particles);
	static std::vector<int> test_kernel(std::vector<Particle>* particles, int watchedId);
	static void test_kernel_integral();
};