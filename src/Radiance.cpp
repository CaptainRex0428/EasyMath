#include "Radiance.h"

#include "EMConst.h"

double EM::RadiansToDegree(double radians)
{
	return (radians * 360) / (2 * PI);
}

float EM::RadiansToDegree(float radians)
{
	return (float)RadiansToDegree((double)radians);
}
