#pragma once

#include "core_math.hpp"



// Scalar equation
constexpr double A = 0.;
constexpr double B = 1.;

inline double R1(double x) {
	return std::cos(std::pow(x, 5) - x + 3 + std::cbrt(2.));
}

inline double R2(double x) {
	return std::atan((cube(x) - 5. * std::sqrt(2) * x - 4) / (std::sqrt(6) * x + std::sqrt(3)));
}

inline double VarF(double x) {
	return R1(x) + R2(x) + 1.8;
}

inline double VarF2(double x) {
	return sqr(VarF(x));
}

inline double VarF3(double x) {
	return cube(VarF(x));
}

inline double AnalythicalDerivative(double x) {
	using namespace std;
	return
		( ((3. * sqr(x) - 5 * sqrt(2)) / (sqrt(6) * x + sqrt(3))) - sqrt(6) * (cube(x) - 5. * sqrt(2) * x - 4) / sqr(sqrt(6) * x + sqrt(3)) )
		/
		(sqr(cube(x) - 5. * sqrt(2) * x - 4) / sqr(sqrt(6) * x + sqrt(3)) + 1)
		-
		(5. * x * x * x * x - 1) * sin(x * x * x * x * x - x + cbrt(2) + 3.);
}

// Analythical derivative

// System of equations
inline double f1(const Vector &X) {
	return sqr(sqr(X(0)) + 2. * sqr(X(1))) - 7. * (sqr(X(0)) - 2. * sqr(X(1)));
}
inline double f2(const Vector &X) {
	return 3. * sqr(X(0)) + 4. * sqr(X(1)) - 14.;
}
inline Vector F(const Vector &X) {
	auto res = Vector{ f1(X), f2(X) };
	return res;
}

// Analythical derivative
inline Matrix AnalythicalJacobian(const Vector &X) {
	return Matrix{
		{-14. * X(0) + 4. * X(0) * (sqr(X(0)) + 2. * sqr(X(1))), 28. * X(1) + 8. * X(1) * (sqr(X(0)) + 2. * sqr(X(1)))},
		{6. * X(0), 8. * X(1)}
	};
}