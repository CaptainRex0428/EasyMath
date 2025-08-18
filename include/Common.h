#pragma once

#include "EasyMathAPI.h"
#include <type_traits>

#define NEARZERO_THRESHOLD 1e-6

namespace EM
{
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