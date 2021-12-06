#pragma once

#include <tuple>

#include "core_math.hpp"


// @return 1 => Solution
// @return 2 => Iterations
// @return 3 => convergence orders
inline std::tuple<double, uint, dvector> bisection_solve(ScalarFunction *f, double left, double right, double precision, uint maxIterations = 50) {
	dvector approx; // vector of aproximations
	dvector orders;

	uint iterations = 0;

	// Get approximatins
	while (true) {
		const double mid = middle(left, right);
		approx.push_back(mid);

		if (sign(f(mid)) == sign(f(left)))
			left = mid; 
		else
			right = mid;

		++iterations;
		if (right - left < 2. * precision || iterations > maxIterations) break;
	}

	// Compute convergence orders
	const double x = middle(left, right);
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
