#include "Vector.h"

#include <utility>
#include <math.h>

#include "Common.h"
#include "Radiance.h"
#include "Degree.h"

bool EM::NearZero(const Vector2& src)
{
	return std::abs(src.Length()) < NEARZERO_THRESHOLD;
}

bool EM::NearZero(const Vector3& src)
{
	return std::abs(src.Length()) < NEARZERO_THRESHOLD;
}

EM::Vector2 EM::DegreeToRadians(const Vector2& degree)
{
	return Vector2{ DegreeToRadians(degree[x]), DegreeToRadians(degree[y]) };
}

EM::Vector3 EM::DegreeToRadians(const Vector3& degree)
{
	return Vector3{ DegreeToRadians(degree[x]), DegreeToRadians(degree[y]), DegreeToRadians(degree[z]) };
}

EM::Vector2 EM::RadiansToDegree(Vector2& radians)
{
	return Vector2{ RadiansToDegree(radians[x]), RadiansToDegree(radians[y]) };
}

EM::Vector3 EM::RadiansToDegree(Vector3& radians)
{
	return Vector3{ RadiansToDegree(radians[x]), RadiansToDegree(radians[y]),RadiansToDegree(radians[z]) };
}
