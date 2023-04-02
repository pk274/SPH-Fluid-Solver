// Paul Kull, 2022
#pragma once

class Parameters {
public:
	static const int WINDOW_WIDTH;
	static const int WINDOW_HEIGHT;
	static const float NEIGHBORHOOD_RADIUS;
	static const int NUM_SUPPOSED_NEIGHBORS;
	static const float GRAVITY;
	static const float GAMMA;
	static const float VISCOSITY;
	static const float H;
	static const float MAX_DENSITY_ERROR;
	static const float OMEGA;
	static const int MAX_SOLVER_ITERATIONS;
	static const float timeStepSize;

	static const float GRAPH_ZOOM;
	static const float GRAPH_SPEED;

	static const int SIMULATION_LENGTH;
	static const bool JUST_RENDER;
	static const bool INTERACTIVE;
	static const int SPEEDUP;
	static const int RENDER_SPEEDUP;
};