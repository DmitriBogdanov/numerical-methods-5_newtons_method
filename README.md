# Numerical methods 5 / Newton's method

Contains implementations of following numerical methods:

* Newton's method for solving nonlinear equations
* Newton's method for solving systems of nonlinear equations

Note that present implementations are intended for academic purposes, as such they are not meant to be used in any sort of high-performance production code.

## Compilation

* Recommended compiler: MSVC v142
* Requires C++17 support

## Dependencies

* Eigen 3.4.0.

## Usage

Place config file of the following format into the same folder as executable:

* Line 1: INPUT_FILEPATH [value without whitespaces]
* Line 2: OUTPUT_FOLDER [value without whitespaces]
* Line 3: PRECISION

Refer to 'config.txt' as an example. Upon execution no furter inputs are required.

## Version history

## License

This project is licensed under the MIT License - see the LICENSE.md file for details