#include "Definitions.hxx"

#ifndef KPM_4_SYSTEMHUBBARD_HXX
#define KPM_4_SYSTEMHUBBARD_HXX


class SystemHubbard3D {
public:
    SystemHubbard3D(int nr, int order, int lx, int ly, int lz, REAL eps, REAL t, REAL I, REAL M_0);
    std::vector<REAL> DOS();
    std::vector<REAL> local_DOS(R3 location);
    std::vector<REAL> SpectralDensity(const std::function<void(vec &, const vec &)> &A);
    std::vector<REAL> local_SpectralDensity(const std::function<void(vec &, const vec &)> &A, R3 loc);
public:
    const int Unit_Cell_size = 2;
    int NR;
    int ORDER;
    int Lx;
    int Ly;
    int Lz;
    int Lsize;
    int SIZE;
    REAL scale;
    REAL T;
    REAL I;
    REAL Chemical_potential;
    std::vector<Neighbor3D> neighborsList;
    std::vector<REAL> M;

    void H(vec &o, const vec &v);
    std::vector<Neighbor3D> GenNeighbors();
    int location(R3 R) const;
};


#endif //KPM_4_SYSTEMHUBBARD_HXX
