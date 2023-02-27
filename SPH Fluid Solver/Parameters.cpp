// Paul Kull, 2022

#include "./Parameters.h"

const int Parameters::WINDOW_WIDTH = 1920;
const int Parameters::WINDOW_HEIGHT = 1080;
const double Parameters::NEIGHBORHOOD_RADIUS = 4;
const int Parameters::NUM_SUPPOSED_NEIGHBORS = 9;
const double Parameters::GRAVITY = 981;
const double Parameters::STIFFNESS = 1000000;	// 1000000
const double Parameters::VISCOSITY = 45;

const float Parameters::H = 2;				// Distance of 2*H is supported by kernel -> H = neigRad / 2
const float Parameters::timeStepSize = 0.001;

const float Parameters::GRAPH_ZOOM = 17;
const float Parameters::GRAPH_SPEED = 2;

const bool Parameters::INTERACTIVE = 1;
const bool Parameters::JUST_RENDER = 0;
const int Parameters::SIMULATION_LENGTH = 6000;
const int Parameters::SPEEDUP = 10;
const int Parameters::RENDER_SPEEDUP = 3;
