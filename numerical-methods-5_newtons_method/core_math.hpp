#pragma once

#include <algorithm>

#include "Eigen/Dense"


using uint = unsigned int;
using sint = int;
using matrix = Eigen::MatrixXd;



constexpr double INF = std::numeric_limits<double>::infinity();



template<typename T>
constexpr bool is_zero(T value) {
	constexpr T EPSILON = static_cast<T>(1e-16);

	return std::abs(value) < EPSILON;
}

// sign(x) = { -1, 0, 1 } signum
template<typename T>
constexpr T sign(T value) {
	return (T(0) < value) - (value < T(0));
}

// sqr(x) = x^2
template<typename T>
constexpr T sqr(T value) {
	return value * value;
}

// cube(x) = x^3
template<typename T>
constexpr T cube(T value) {
	return value * value * value;
}
