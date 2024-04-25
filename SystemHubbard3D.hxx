#include "Definitions.hxx"

#ifndef KPM_4_SYSTEMHUBBARD_HXX
#define KPM_4_SYSTEMHUBBARD_HXX


class SystemHubbard3D {
public:
    SystemHubbard3D(int nr, int order, int lx, int ly, int lz, REAL eps, REAL t, REAL I, REAL M_0);
    std::vector<REAL> DOS() const;
    std::vector<REAL> local_DOS(R3 location)const;
    std::vector<REAL> SpectralDensity(const std::function<void(vec &, const vec &)> &A) const;
    std::vector<std::vector<REAL>> SpectralDensityMulti(const h_func &A) const;
    std::vector<REAL> LocalSpectralDensity(const std::function<void(vec &, const vec &)> &A, R3 loc) const;
    std::vector<std::vector<REAL>> LocalSpectralDensityMulti(const h_func &A, R3 loc) const;

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

    void H(vec &o, const vec &v) const;
    std::vector<Neighbor3D> GenNeighbors();
    int location(R3 R) const;
    int lat_location(R3 R) const;

};


#endif //KPM_4_SYSTEMHUBBARD_HXX
