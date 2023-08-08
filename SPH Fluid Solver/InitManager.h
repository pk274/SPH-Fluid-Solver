// Paul Kull, 2023
#pragma once

#include <vector>
#include <SFML/Graphics.hpp>
#include <fstream>

#include "./Functions.h"
#include "./Particle.h"
#include "./SolidParticle.h"
#include "./FluidParticle.h"
#include "./HashManager.h"
#include "./TestManager.h"
#include "./CompactHashManager.h"
#include "./Simulation.h"
#include "./Parameters.h"

enum SimulationPreset {
	Alone = 0,
	SmallBox = 1,
	StuffedBox = 2,
	SingleParticle = 3,
	RotatedBox = 4,
	FewParticles = 5,
	GiantBox = 6,
	FourLayers = 7,
	ManyLayers = 8,
	WideLayers = 9,
	BreakingDam = 10,
	BigBreakingDam = 11,
	FixedBreakingDam = 12,
	SmallBreakingDam = 13,
	TallBreakingDam = 14,
	Cup = 20,
	Complex = 30,
	Complex2 = 31,
	Osmosis = 40,
	SideSpawn = 50,
	Rain = 51,
	Rain2 = 52,
	Fountain = 60,
	MovingWall = 70,
	WaveGenerator = 71,
	Island = 75,
	CubeDrop = 80,
	CompressionState = 90,
};


class InitManager {
public:

	Simulation* _sim;


	InitManager(Simulation* sim);

	void init_simulation(SimulationPreset preset);

	void init_empty_simulation();
	void init_stuffed_box_simulation(int size = 50, int zoom = 5);
	void init_single_particle_simulation(int size = 50, int zoom = 5);
	void init_rotated_box_simulation(int size, int zoom, int rotation);
	void init_random_particles_simulation(int size, int zoom, int numParticles);
	void init_breaking_dam_simulation(int size = 50, int zoom = 5, bool adaptive = true, int numFluidParticles = 10000);
	void init_tall_breaking_dam_simulation(int size = 50, int zoom = 5);
	void init_layer_simulation(int size, int zoom, int layers, bool xOffset = true, bool yOffset = true);
	void init_wide_layer_simulation(int size, int zoom, int layers, bool xOffset = true, bool yOffset = true);
	void init_cup_simulation(int size, int zoom);
	void init_complex_simulation(int size, int zoom);
	void init_complex_2_simulation(int size, int zoom);
	void init_osmosis_simulation(int size, int zoom);
	void init_side_spawn_simulation(int size, int zoom);
	void init_rain_simulation(int size, int zoom);
	void init_rain2_simulation(int size, int zoom);
	void init_fountain_simulation();
	void init_moving_wall_simulation(int size, int zoom);
	void init_wave_generator_simulation(int size, int zoom);
	void init_island_simulation(int size, int zoom);
	void init_cube_drop_simulation(int size, int zoom);
	void init_compressed_state_simulation(int size, int zoom);
};
