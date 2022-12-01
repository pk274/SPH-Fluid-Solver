// Paul Kull, 2022
#pragma once

class Parameters {
public:
	static const int WINDOW_WIDTH;
	static const int WINDOW_HEIGHT;
	static const double NEIGHBORHOOD_RADIUS;
	static const int NUM_SUPPOSED_NEIGHBORS;
	static const double GRAVITY;
	static const double STIFFNESS;
	static const double VISCOSITY;
	static const float H;
	static const float timeStepSize;

	static const float GRAPH_ZOOM;
	static const float GRAPH_SPEED;

	static const int SIMULATION_LENGTH;
	static const bool JUST_RENDER;
	static const bool INTERACTIVE;
	static const int SPEEDUP;
};