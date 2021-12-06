#pragma once

#include <tuple>
#include <vector>

#include "core_math.hpp"



// @return 1 => Solution
// @return 2 => Iterations
// @return 3 => convergence orders
inline std::tuple<double, uint, dvector> secant_solve(ScalarFunction *f, double A, double B, double precision, uint maxIterations = 50) {
	dvector approx; // vector of aproximations
	dvector orders;

	uint iterations = 0;

	double x0 = A;
	double x = x0;
	approx.push_back(x0);

	while (true) {
		x = x0 - f(x0) * (B - x0) / (f(B) - f(x0));
		approx.push_back(x);

		++iterations;
		if (std::abs(x - x0) < precision || iterations > maxIterations) break;

		x0 = x;
	}

	// Compute convergence orders
	const auto N = approx.size();
	orders.reserve(N > 2 ? N - 2 : 0);

	for (size_t k = 1; k < N - 2; ++k) {
		const double err_ = std::abs(approx[k - 1] - x);
		const double _err_ = std::abs(approx[k] - x);
		const double _err = std::abs(approx[k + 1] - x);
		const double p = std::log(_err / _err_) / std::log(_err_ / err_);
		orders.push_back(p);
	}

	return { x, iterations, orders };
}



// @return 1 => Solution
// @return 2 => Iterations
// @return 3 => convergence orders
inline std::tuple<double, uint, dvector> newton_solve(ScalarFunction *f, double x0, double precision, bool dIsNumeric = true, uint maxIterations = 50) {
	dvector approx; // vector of aproximations
	dvector orders;

	uint iterations = 0;

	double x = x0;
	approx.push_back(x0);

	while (true) {
		const auto Derivative = dIsNumeric ? derivative(f, x0) : AnalythicalDerivative(x0);

		x = x0 - f(x0) / Derivative;
		approx.push_back(x);

		++iterations;
		if (std::abs(x - x0) < precision || iterations > maxIterations) break;

		x0 = x;
	}

	// Compute convergence orders
	const auto N = approx.size();
	orders.reserve(N > 2 ? N - 2 : 0);

	for (size_t k = 1; k < N - 2; ++k) {
		const double err_ = std::abs(approx[k - 1] - x);
		const double _err_ = std::abs(approx[k] - x);
		const double _err = std::abs(approx[k + 1] - x);
		const double p = std::log(_err / _err_) / std::log(_err_ / err_);
		orders.push_back(p);
	}

	return { x, iterations, orders };
}


// @return 1 => Solution
// @return 2 => Iterations
// @return 3 => convergence orders
inline std::tuple<Vector, uint, dvector> newton_solve_system(VectorFunction *F, Vector X0, double precision, bool dIsNumeric = true, uint maxIterations = 50) {
	vvector approx; // vector of aproximations
	dvector orders;

	uint iterations = 0;
	
	Vector X = X0;
	approx.push_back(X0);

	// Fill approximations
	while (true) {
		const auto Jacobian = dIsNumeric ? jacobian(F, X0) : AnalythicalJacobian(X0);

		X = X0 - Jacobian.inverse() * F(X0);
		approx.push_back(X);

		++iterations;
		if ((X - X0).norm() < precision || iterations > maxIterations) break;

		X0 = X;
	}

	// Compute convergence orders
	const auto N = approx.size();
	orders.reserve(N > 2 ? N - 2 : 0);

	for (size_t k = 1; k < N - 2; ++k) {
		const double err_ = (approx[k - 1] - X).norm();
		const double _err_ = (approx[k] - X).norm();
		const double _err = (approx[k + 1] - X).norm();
		const double p = std::log(_err / _err_) / std::log(_err_ / err_);
		orders.push_back(p);
	}

	return { X, iterations, orders };
}