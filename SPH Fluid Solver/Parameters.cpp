// Paul Kull, 2022

#include <SFML/Graphics.hpp>
#include "./Parameters.h"

const bool Parameters::SOLVE_PPE = 1;
const bool Parameters::S_DI = 0;
const bool Parameters::S_VD = 0;
const bool Parameters::S_VD_DI = 0;
const bool Parameters::S_VD_DI_SIMPLE = 1;

const int Parameters::WINDOW_WIDTH = 1920;
const int Parameters::WINDOW_HEIGHT = 1080;
const float Parameters::NEIGHBORHOOD_RADIUS = 4;
const int Parameters::NUM_SUPPOSED_NEIGHBORS = 9;
const float Parameters::GRAVITY = 981; // 981;
const float Parameters::GAMMA = 0.7;
const float Parameters::BOUNDARY_VISCOSITY = 3.f;		// 0.3

const float Parameters::H = 2;				// Distance of 2*H is supported by kernel -> H = neigRad / 2
const float Parameters::MAX_DENSITY_ERROR = 0.001;	// vd 0.00001
const float Parameters::OMEGA = 0.5;
const int Parameters::MAX_SOLVER_ITERATIONS = 100;
const int Parameters::MIN_SOLVER_ITERATIONS = 2;
const float Parameters::TIME_STEP = 0.001;	//VD 0.002

const float Parameters::GRAPH_ZOOM = 17;
const float Parameters::GRAPH_SPEED = 2;

const bool Parameters::INTERACTIVE = 1;
const bool Parameters::JUST_RENDER = 0;
const float Parameters::SIMULATION_LENGTH = 0.2f;

const bool Parameters::COLOR_CODE_SPEED = 1;
const bool Parameters::COLOR_CODE_PRESSURE = 0;
const bool Parameters::COLOR_CODE_DENSITY = 0;

const bool Parameters::ADAPTIVE_TIME_STEP = 0;
const float Parameters::MAX_TIME_STEP = 0.006;
const float Parameters::INITIALIZATION_PHASE = 0.1f;
const float Parameters::CFL_NUMBER = 0.7;
const int Parameters::DESIRED_ITERATIONS = 17;
const float Parameters::ITERATIONS_ENFORCEMENT = 0.6;	// (0, 1]
const double Parameters::TIME_OFFSET = 0.0000001;
const float Parameters::SLOW_DOWN = 1;
const float Parameters::SMOOTHING = 0.f;

const float Parameters::PRESSURE_CODE_ROUGHNESS = 1.5;
const float Parameters::VELOCITY_CODE_RANGE = 1.5;
const float Parameters::DENSITY_CODE_INTENSITY = 100;

const bool Parameters::DOCUMENT_AVG_DENSITY = 1;
const bool Parameters::DOCUMENT_ITERATIONS = 1;
const bool Parameters::DOCUMENT_TIME = 1;
const bool Parameters::DOCUMENT_ESTIMATED_DENSITY = 1;
const bool Parameters::WRITE_SCREEN_IMAGES = 0;

const float Parameters::EOS_STIFFNESS = 800000;	// 1000000
const float Parameters::EOS_VISCOSITY = 35;

const sf::Color Parameters::BACKGROUND_COLOR = sf::Color::Black;	//230
const sf::Color Parameters::SOLID_COLOR = sf::Color::White;		//150
const bool Parameters::PRETTY_MODE = 0;

//const sf::Color Parameters::BACKGROUND_COLOR = sf::Color::Color(230, 230, 230);	//230
//const sf::Color Parameters::SOLID_COLOR = sf::Color::Color(150, 150, 150);		//150
