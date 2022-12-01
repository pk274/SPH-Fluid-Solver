// Paul Kull, 2022

#include "./Parameters.h"

const int Parameters::WINDOW_WIDTH = 1100;
const int Parameters::WINDOW_HEIGHT = 800;
const double Parameters::NEIGHBORHOOD_RADIUS = 4;
const int Parameters::NUM_SUPPOSED_NEIGHBORS = 7;
const double Parameters::GRAVITY = 980;
const double Parameters::STIFFNESS = 1000000;	
const double Parameters::VISCOSITY = 35;	// 5000 so funny
const float Parameters::H = 2;				// Distance of 2*H is supported by kernel -> H = neigRad / 2
const float Parameters::timeStepSize = 0.00005;		// 0.00007

const float Parameters::GRAPH_ZOOM = 6;
const float Parameters::GRAPH_SPEED = 2;

const bool Parameters::INTERACTIVE = 0;
const bool Parameters::JUST_RENDER = 0;
const int Parameters::SIMULATION_LENGTH = 50000;
const int Parameters::SPEEDUP = 10;
