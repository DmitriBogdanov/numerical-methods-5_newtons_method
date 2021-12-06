#include <iostream>
#include <iomanip>
#include <fstream>

#include "equations.hpp"
#include "bisection.hpp"
#include "newton.hpp"



const std::string INPATH = "[input]/";
const std::string OUTPATH = "[output]/";
constexpr double PRECISION = 1e-6;

constexpr size_t GRID_SIZE = 100;
constexpr double GRID_X1_MIN = -6;
constexpr double GRID_X1_MAX = 6;
constexpr double GRID_X2_MIN = -6;
constexpr double GRID_X2_MAX = 6;
const double STEP_X1 = (GRID_X1_MAX - GRID_X1_MIN) / (GRID_SIZE - 0);
const double STEP_X2 = (GRID_X2_MAX - GRID_X2_MIN) / (GRID_SIZE - 0);
constexpr uint NEWTON_MAX_ITERATIONS = 30;

constexpr std::streamsize SETW_SIZE = 19;

void run_bisection() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Bysection running...";

	const auto [solution, iterations, orders] = bisection_solve(VarF, A, B, PRECISION);;

	std::cout
		<< STARTER << "Method     -> Bisection"
		<< STARTER << "Precision  -> " << PRECISION
		<< STARTER << "Interval   -> [" << A << ", " << B << "]"
		<< STARTER << "Solution   -> " << solution
		<< STARTER << "Iterations -> " << iterations;

	// Save to files
	std::ofstream fileAB(OUTPATH + "[bisection](AB).txt");
	fileAB << A << ' ' << B;

	std::ofstream filePrecision(OUTPATH + "[bisection](precision).txt");
	filePrecision << PRECISION;

	std::ofstream fileSolution(OUTPATH + "[bisection](solution).txt");
	fileSolution << solution;

	std::ofstream fileIterations(OUTPATH + "[bisection](iterations).txt");
	fileIterations << iterations;

	std::ofstream fileOrders(OUTPATH + "[bisection](orders).txt");
	for (const auto &el : orders) fileOrders << el << '\n';
}

void run_secant() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Secant running...";

	const auto [solution, iterations, orders] = secant_solve(VarF, A, B, PRECISION);;

	std::cout
		<< STARTER << "Method     -> Secant"
		<< STARTER << "Precision  -> " << PRECISION
		<< STARTER << "Interval   -> [" << A << ", " << B << "]"
		<< STARTER << "Solution   -> " << solution
		<< STARTER << "Iterations -> " << iterations;

	// Save to files
	std::ofstream fileAB(OUTPATH + "[secant](AB).txt");
	fileAB << A << ' ' << B;

	std::ofstream filePrecision(OUTPATH + "[secant](precision).txt");
	filePrecision << PRECISION;

	std::ofstream fileSolution(OUTPATH + "[secant](solution).txt");
	fileSolution << solution;

	std::ofstream fileIterations(OUTPATH + "[secant](iterations).txt");
	fileIterations << iterations;

	std::ofstream fileOrders(OUTPATH + "[secant](orders).txt");
	for (const auto &el : orders) fileOrders << el << '\n';
}

void run_scalar_newton(bool dIsNumeric) {
	std::cout
		<< SEPARATOR
		<< STARTER << "Newton's (scalar) running...";

	const auto x0 = 0.5;
	const auto [solution, iterations, orders] = newton_solve(VarF3, x0, PRECISION, dIsNumeric);

	std::cout
		<< STARTER << "Method     -> Newton's method (scalar function)"
		<< STARTER << "Precision  -> " << PRECISION
		<< STARTER << "Initial X  -> " << x0
		<< STARTER << "dIsNumeric -> " << dIsNumeric
		<< STARTER << "Solution   -> " << solution
		<< STARTER << "Iterations -> " << iterations;

	// Save to files
	std::ofstream fileX0(OUTPATH + "[newton](X0).txt");
	fileX0 << x0;

	std::ofstream filePrecision(OUTPATH + "[newton](precision).txt");
	filePrecision << PRECISION;

	std::ofstream fileDerivative(OUTPATH + "[newton](derivative).txt");
	fileDerivative << dIsNumeric;

	std::ofstream fileSolution(OUTPATH + "[newton](solution).txt");
	fileSolution << solution;

	std::ofstream fileIterations(OUTPATH + "[newton](iterations).txt");
	fileIterations << iterations;

	std::ofstream fileOrders(OUTPATH + "[newton](orders).txt");
	for (const auto &el : orders) fileOrders << el << '\n';
}

void run_vector_newton(bool dIsNumeric) {
	std::cout
		<< SEPARATOR
		<< STARTER << "Newton's (vector) running...";

	const auto X0 = Vector{ 3., 1. };
	//const auto X0 = Vector::Zero();
	const auto [solution, iterations, orders] = newton_solve_system(F, X0, PRECISION, dIsNumeric);

	std::cout
		<< STARTER << "Method     -> Newton's method (vector function)"
		<< STARTER << "Precision  -> " << PRECISION
		<< STARTER << "Initial X  -> " << X0.format(INLINE)
		<< STARTER << "dIsNumeric -> " << dIsNumeric
		<< STARTER << "Solution   -> " << solution.format(INLINE)
		<< STARTER << "Iterations -> " << iterations;

	// Save to files
	std::ofstream fileX0(OUTPATH + "[vector_newton](X0).txt");
	fileX0 << X0(0) << ' ' << X0(1);

	std::ofstream filePrecision(OUTPATH + "[vector_newton](precision).txt");
	filePrecision << PRECISION;

	std::ofstream fileDerivative(OUTPATH + "[vector_newton](derivative).txt");
	fileDerivative << dIsNumeric;

	std::ofstream fileSolution(OUTPATH + "[vector_newton](solution).txt");
	fileSolution << solution(0) << ' ' << solution(1);

	std::ofstream fileIterations(OUTPATH + "[vector_newton](iterations).txt");
	fileIterations << iterations;

	std::ofstream fileOrders(OUTPATH + "[vector_newton](orders).txt");
	for (const auto &el : orders) fileOrders << el << '\n';
}

Vector gridPos(size_t i, size_t j) {
	return Vector{ GRID_X1_MIN + i * STEP_X1, GRID_X2_MIN + j * STEP_X2 };
}

void run_newton_convergence_area_test() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Newton's convergence area estimate running...";

	// Evaluate at grid, save number of iteration at each point in a single long list
	std::ofstream filePoints(OUTPATH + "[points].txt");
	std::ofstream fileInfo(OUTPATH + "[info].txt");

	constexpr double PROGRESSBAR_STEP = 0.1;
	
	double progressbarCounter = 0;
	uint progressbarPercentage = 0;

	for (size_t i = 0; i <= GRID_SIZE; ++i) {
		for (size_t j = 0; j <= GRID_SIZE; ++j) {
			const auto point = gridPos(i, j);
			const auto iterations = std::get<1>(newton_solve_system(F, gridPos(i, j), PRECISION, NEWTON_MAX_ITERATIONS));

			filePoints
				<< std::setw(SETW_SIZE) << point(0)
				<< std::setw(SETW_SIZE) << point(1)
				<< std::setw(SETW_SIZE) << iterations
				<< '\n';
		}

		progressbarCounter += 1. / (GRID_SIZE + 1);

		if (progressbarCounter > PROGRESSBAR_STEP) {
			progressbarCounter -= PROGRESSBAR_STEP;
			progressbarPercentage += static_cast<uint>(100 * PROGRESSBAR_STEP);
			std::cout << STARTER << progressbarPercentage << "% done...";
		}
	}

	fileInfo
		<< std::setw(SETW_SIZE) << GRID_X1_MIN << '\n'
		<< std::setw(SETW_SIZE) << GRID_X1_MAX << '\n'
		<< std::setw(SETW_SIZE) << GRID_X2_MIN << '\n'
		<< std::setw(SETW_SIZE) << GRID_X2_MAX << '\n'
		<< std::setw(SETW_SIZE) << NEWTON_MAX_ITERATIONS << '\n';

	std::cout
		<< STARTER << "Evaluation finished. Results saved to corresponding files.";
}

int main(int argc, char* argv[]) {
	// Try methods one by one
	run_bisection();
	run_secant();

	constexpr bool dIsNumeric = true;
	run_scalar_newton(dIsNumeric);
	run_vector_newton(dIsNumeric);

	// Find newton's method convergence area
	run_newton_convergence_area_test();

	std::cout << "\n\n>>> Execution finished successfully.\n";

	return 0;
}