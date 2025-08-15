#include "Common.h"

#include <cmath>

int EM::Clamp(int src, int min, int max)
{
	if (src <= min)
	{
		return min;
	}

	if (src >= max)
	{
		return max;
	}

	return src;
}

float EM::Clamp(float src, float min, float max)
{
	if (src <= min)
	{
		return min;
	}

	if (src >= max)
	{
		return max;
	}

	return src;
}

uint8_t EM::Clamp(uint8_t src, uint8_t min, uint8_t max)
{
	if (src <= min)
	{
		return min;
	}

	if (src >= max)
	{
		return max;
	}

	return src;
}

double& EM::limit_min(double& src, double min)
{
	if (src < min)
	{
		src = min;
	}

	return src;
}

float& EM::limit_min(float& src, float min)
{
	if (src < min)
	{
		src = min;
	}

	return src;
}

int& EM::limit_min(int& src, int min)
{
	if (src < min)
	{
		src = min;
	}

	return src;
}

size_t& EM::limit_min(size_t& src, size_t min)
{
	if (src < min)
	{
		src = min;
	}

	return src;
}

double& EM::limit_max(double& src, double max)
{
	if (src > max)
	{
		src = max;
	}

	return src;
}

float& EM::limit_max(float& src, float max)
{
	if (src > max)
	{
		src = max;
	}

	return src;
}

int& EM::limit_max(int& src, int max)
{
	if (src > max)
	{
		src = max;
	}

	return src;
}

size_t& EM::limit_max(size_t& src, size_t max)
{
	if (src > max)
	{
		src = max;
	}

	return src;
}

bool EM::NearZero(double& src)
{
	return NearZero(std::move(src));
}

bool EM::NearZero(double&& src)
{
	return std::abs(src) < NEARZERO_THRESHOLD;
}

bool EM::NearZero(float& src)
{
	return NearZero(std::move(src));
}

bool EM::NearZero(float&& src)
{
	return NearZero(double(src));
}