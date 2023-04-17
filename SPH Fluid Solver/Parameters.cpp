// Paul Kull, 2022

#include "./Parameters.h"

const int Parameters::WINDOW_WIDTH = 1920;
const int Parameters::WINDOW_HEIGHT = 1080;
const float Parameters::NEIGHBORHOOD_RADIUS = 4;
const int Parameters::NUM_SUPPOSED_NEIGHBORS = 9;
const float Parameters::GRAVITY = 981; // 981;
const float Parameters::GAMMA = 0.7;	// 1000000
const float Parameters::BOUNDARY_VISCOSITY = 3;

const float Parameters::H = 2;				// Distance of 2*H is supported by kernel -> H = neigRad / 2
const float Parameters::MAX_DENSITY_ERROR = 0.01;
const float Parameters::OMEGA = 0.5;
const int Parameters::MAX_SOLVER_ITERATIONS = 10;
const float Parameters::TIME_STEP = 0.003;

const float Parameters::GRAPH_ZOOM = 17;
const float Parameters::GRAPH_SPEED = 2;

const bool Parameters::INTERACTIVE = 1;
const bool Parameters::JUST_RENDER = 1;
const int Parameters::SIMULATION_LENGTH = 10;

const bool Parameters::COLOR_CODE_SPEED = 1;
const bool Parameters::COLOR_CODE_PRESSURE = 0;
const bool Parameters::ADAPTIVE_TIME_STEP = 1;
const float Parameters::MAX_TIME_STEP = 0.005;
const float Parameters::CFL_NUMBER = 0.7;
const double Parameters::TIME_OFFSET = 0.000001;
const float Parameters::SLOW_DOWN = 1;

const float Parameters::PRESSURE_CODE_ROUGHNESS = 2;