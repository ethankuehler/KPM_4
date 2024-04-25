#include <iostream>
#include <array>
#include <complex>
#include <random>
#include <fstream>
#include <iomanip>
#include </usr/local/include/eigen3/Eigen/Core>

#ifndef KPM_4_DEFINITIONS_HXX
#define KPM_4_DEFINITIONS_HXX

typedef double REAL;

using namespace std::complex_literals;

//3d lattice vector.
typedef std::array<int, 3> R3;
//2d lattice vector.
typedef std::array<int, 2> R2;
//3D Neighbor +x, -x, +y, -y, +z, -z
typedef std::array<int, 7> Neighbor3D;
//2D Neighbor +x, -x, +y, -y
typedef std::array<int, 5> Neighbor2D;
//Vector of arbitrary size, complex double.
typedef Eigen::Matrix<std::complex<REAL>, -1, 1> vec;
//Lattice site;
typedef Eigen::Matrix<std::complex<REAL>, 2, 1> site2D;
//2D matrix
typedef Eigen::Matrix<std::complex<REAL>, 2, 2> mat2D;
//complex number
typedef std::complex<REAL> COMPLEX;
//1i but setup for different types of REAL
const std::complex<REAL> ONE{0, 1};
//H function
typedef std::vector<std::function<void(vec &, const vec &)>> h_func;

#endif //KPM_4_DEFINITIONS_HXX
