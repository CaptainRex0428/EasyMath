#pragma once

#include "EasyMathAPI.h"
#include <type_traits>

#define NEARZERO_THRESHOLD 0.0001

namespace EM
{
	EASYMATH_API int Clamp(int src, int min = 0, int max = 1);
	EASYMATH_API float Clamp(float src, float min = 0, float max = 1);
	EASYMATH_API uint8_t Clamp(uint8_t src, uint8_t min = 0, uint8_t max = 1);

	EASYMATH_API double& limit_min(double& src, double min);
	EASYMATH_API float& limit_min(float& src, float min);
	EASYMATH_API int& limit_min(int& src, int min);
	EASYMATH_API size_t& limit_min(size_t& src, size_t min);

	EASYMATH_API double& limit_max(double& src, double max);
	EASYMATH_API float& limit_max(float& src, float max);
	EASYMATH_API int& limit_max(int& src, int max);
	EASYMATH_API size_t& limit_max(size_t& src, size_t max);

	EASYMATH_API bool NearZero(double& src);
	EASYMATH_API bool NearZero(double&& src);
	EASYMATH_API bool NearZero(float& src);
	EASYMATH_API bool NearZero(float&& src);

	template<typename T, typename... Args>
	struct are_comparable :
		std::conjunction<std::is_convertible<decltype(std::declval<T>() < std::declval<Args>()), bool>...> {};

	template<typename First, typename... Rest,
		typename = std::enable_if_t<are_comparable<First, Rest...>::value>>
		auto Max(const First& first, const Rest&... rest) {
		if constexpr (sizeof...(rest) == 0) {
			return first;
		}
		else {
			auto max_rest = Max(rest...);
			return first > max_rest ? first : max_rest;
		}
	}

	template<typename First, typename... Rest,
		typename = std::enable_if_t<are_comparable<First, Rest...>::value>>
		auto Min(const First& first, const Rest&... rest) {
		if constexpr (sizeof...(rest) == 0) {
			return first;
		}
		else {
			auto min_rest = Min(rest...);
			return first < min_rest ? first : min_rest;
		}
	}


}