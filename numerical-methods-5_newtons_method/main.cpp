#include <iostream>

#include "core_math.hpp"

constexpr double A = 0.;
constexpr double B = 1.;

double R1(double x) {
	return std::cos( std::pow(x, 5) - x + 3 + std::cbrt(2.) );
}

double R2(double x) {
	return std::atan( (cube(x) + 3. * std::sqrt(2) * x - 4) / (std::sqrt(6) * x + std::sqrt(3)) );
}

double VarF(double x) {
	return R1(x) + R2(x) + 1.8;
}

double f1(double x1, double x2) {
	return std::exp(x1) - x2 + 4.;
}
double f2(double x1, double x2) {
	return sqr(x1) - 2. * x2 - 3.;
}


int main(int argc, char* argv[]) {


	return 0;
}