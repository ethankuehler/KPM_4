#include "Definitions.hxx"

#ifndef KPM_4_SYSTEMHUBBARD2D_HXX
#define KPM_4_SYSTEMHUBBARD2D_HXX


class SystemHubbard2D {
public:
    SystemHubbard2D(int nr, int order, int lx, int ly, REAL eps, REAL t, REAL U, REAL M_0);
    std::vector<REAL> DOS() const;
    std::vector<REAL> local_DOS(R2 loc) const;
    std::vector<std::vector<REAL>> SpectralDensity(std::vector<std::function<void(vec &, const vec &)>>& A) const;
    std::vector<REAL> SpectralDensity(const std::function<void(vec &, const vec &)> &A) const;
    std::vector<REAL> local_SpectralDensity(const std::function<void(vec &, const vec &)> &A, R2 loc) const;

    const int Unit_Cell_size = 2;
    int NR;
    int ORDER;
    int Lx;
    int Ly;
    int Lsize;
    int SIZE;
    REAL scale;
    REAL T;
    REAL I;
    REAL Chemical_potential;
    std::vector<Neighbor2D> neighborsList;
    std::vector<REAL> M;

    void H(vec &o, const vec &v) const;

    std::vector<Neighbor2D> GenNeighbors() const;

    int location(R2 R) const;
};



#endif //KPM_4_SYSTEMHUBBARD2D_HXX
