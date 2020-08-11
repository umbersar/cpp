#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include<numbers>


// Usings

using std::shared_ptr;
using std::make_shared;
using std::sqrt;

// Constants

const double infinity = std::numeric_limits<double>::infinity();
//const double pi = 3.1415926535897932385;
const double pi = std::numbers::pi;

// Utility Functions

inline double degrees_to_radians(double degrees) {
	return degrees * pi / 180.0;
}

//inline double random_double(double lowerLimit = 0, double upperLimit = 1) {
//	std::uniform_real_distribution<double> dist(lowerLimit, upperLimit);  //(min, max)
//	//Mersenne Twister: Good quality random number generator
//	std::mt19937 rng;
//	//Initialize with non-deterministic seeds
//	rng.seed(std::random_device{}());
//	return dist(rng);
//}

inline double clamp(double x, double min, double max) {
	//if (x < min) return min;
	//if (x > max) return max;
	//return x;
	
	return (x < min) ? min : ((x > max) ? max : x);
}

// Common Headers

#include "ray.h"
#include "vec3.h"

#endif