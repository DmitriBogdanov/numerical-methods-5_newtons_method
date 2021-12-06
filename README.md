# Numerical methods 5 / Newton's method

Contains implementations of following numerical methods:

* Bisection method
* Secant method
* Newton's method
* Newton's method for systems

Note that present implementations are intended for academic purposes, as such they are not meant to be used in any sort of high-performance production code.

## Compilation

* Recommended compiler: MSVC v142
* Requires C++17 support

## Dependencies

* Eigen 3.4.0.

## Usage

Precision, max iterations and other parameters are entered as consts before main(). 

## Version history

* 00.02
    * Implemented Newton's method
    * Implemented Newton's method for systems
    * Implemented secant method
    * Implemented bisection method
    * Both Newton methods can be toggled between analythical and numeric derivatives
    * Impremented trackers for convergence order in each method

* 00.01
    * Created core math methods

## License

This project is licensed under the MIT License - see the LICENSE.md file for details