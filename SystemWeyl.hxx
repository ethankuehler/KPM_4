#include "Definitions.hxx"

#ifndef KPM_4_SYSTEMWEYL_HXX
#define KPM_4_SYSTEMWEYL_HXX


class SystemWeyl {
public:
    SystemWeyl(int nr, int order, int lx, int ly, int lz, REAL eps, REAL t, REAL tp, REAL m);
    std::vector<REAL> DOS();
    std::vector<REAL> local_DOS(R3 loc);
    std::vector<REAL> SpectralDensity(const std::function<void(vec &, const vec &)> &A);
    std::vector<REAL> local_SpectralDensity(const std::function<void(vec &, const vec &)> &A);
private:
    const int Unit_Cell_size = 2;
    int NR;
    int ORDER;
    int Lx;
    int Ly;
    int Lz;
    int Lsize;
    int SIZE;
    REAL scale;
    REAL t; // already take into account t/2
    REAL tp;
    REAL m;
    std::vector<Neighbor3D> neighborsList;

    void H(vec &o, const vec &v);
    std::vector<Neighbor3D> GenNeighbors();
    int location(R3 R) const;
};


#endif //KPM_4_SYSTEMWEYL_HXX
