#include "Degree.h"

#include "EMConst.h"

double EM::DegreeToRadians(double degree)
{
	return degree * 2 * PI / 360;
}

float EM::DegreeToRadians(float degree)
{
	return (float)DegreeToRadians((double)degree);
}