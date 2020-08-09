#pragma once
#include"vec3.h"

class ray
{
private:
	point3 orig;
	vec3 dir;
public:
	ray();
	ray(const point3& origin, const vec3& direction);

	point3 origin() const;
	point3 direction() const;
	point3 at(double t);
};

