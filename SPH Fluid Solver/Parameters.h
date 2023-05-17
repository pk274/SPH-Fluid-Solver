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
	static const float BOUNDARY_VISCOSITY;
	static const float H;
	static const float MAX_DENSITY_ERROR;
	static const float OMEGA;
	static const int MAX_SOLVER_ITERATIONS;
	static const int MIN_SOLVER_ITERATIONS;
	static const float TIME_STEP;

	static const float GRAPH_ZOOM;
	static const float GRAPH_SPEED;

	static const int SIMULATION_LENGTH;
	static const bool JUST_RENDER;
	static const bool INTERACTIVE;

	static const bool COLOR_CODE_SPEED;
	static const bool COLOR_CODE_PRESSURE;
	static const bool COLOR_CODE_DENSITY;
	static const bool ADAPTIVE_TIME_STEP;
	static const float MAX_TIME_STEP;
	static const float INITIALIZATION_PHASE;
	static const float CFL_NUMBER;
	static const double TIME_OFFSET;
	static const float SLOW_DOWN;

	static const float PRESSURE_CODE_ROUGHNESS;
	static const float VELOCITY_CODE_RANGE;
	static const float DENSITY_CODE_INTENSITY;

	static const bool DOCUMENT_ITERATIONS;
	static const bool DOCUMENT_TIME;
	static const bool DOCUMENT_AVG_DENSITY;
	static const bool DOCUMENT_ESTIMATED_DENSITY;
	static const bool WRITE_SCREEN_IMAGES;

	static const float EOS_STIFFNESS;
	static const float EOS_VISCOSITY;

	static const sf::Color BACKGROUND_COLOR;
	static const sf::Color SOLID_COLOR;
};