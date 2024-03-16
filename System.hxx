#include <iostream>
#include <array>
#include <complex>
#include <random>
#include <fstream>
#include <iomanip>
#include </usr/local/include/eigen3/Eigen/Core>
#include "Definitions.hxx"

#ifndef KPM_4_SYSTEM_HXX
#define KPM_4_SYSTEM_HXX


//its +x and -x
typedef std::array<int, 3> Neighbor1D;


struct TightBinding1D {
    TightBinding1D(REAL t, int nr, int order, int lx, REAL eps, REAL s);
    int Unit_Cell_size = 2;
    int NR;
    int ORDER;
    int Lx;
    int SIZE;
    REAL t;
    REAL scale;
    std::vector<Neighbor1D> neighbors;
    std::vector<REAL> DOS();
    std::vector<REAL> Spectral_density(const std::function<void(vec&, const vec&)>& A);

    void gen_neighbors();
    void H(vec &out, const vec &in);
};



#endif //KPM_4_SYSTEM_HXX
