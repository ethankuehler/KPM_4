//
// Created by bruv island on 3/2/24.
//

#include "SystemWeyl.hxx"

SystemWeyl::SystemWeyl(int nr, int order, int lx, int ly, int lz, REAL eps, REAL t, REAL tp, REAL m, REAL i)
: NR(nr), ORDER(order), Lx(lx), Ly(ly), Lz(lz), t(t), tp(tp), m(m), I(i) {
    Lsize = Lx * Ly * Lz;
    SIZE = Lsize*Unit_Cell_size;
    scale = (9+I)/(2 - eps);
    this->t = t/(2*scale);
    this->tp = tp/(2*scale);
    this->m = m/scale;
    this->I = I/scale;
    neighborsList = GenNeighbors();
    M = std::vector<REAL>(Lsize, 0);
}

std::vector<REAL> SystemWeyl::DOS() const {
    std::vector<REAL> mu = std::vector<REAL>(ORDER, 0);

    //random number engine.
    std::random_device rd;  // obtain a random seed from the OS
    std::mt19937 eng(rd());  // seed the generator
    std::uniform_real_distribution<REAL> distr(0, 2*M_PI);

#pragma omp parallel for
    for (int i = 0; i < NR; i++) {
        vec v(SIZE);
        for (auto &j: v) {
            j = std::exp(distr(eng)*ONE);
        }
        v.normalize();
        vec T0 = v;
        vec T1(SIZE);
        H(T1, v);
        vec T2(SIZE);

        std::complex<REAL> weight = 0;
        for (int j = 0; j < ORDER; j++) {
            weight = v.adjoint()*T0;
            mu[j] += weight.real();
            H(T2, 2*T1);
            T2 -= T0;
            T0 = T1;
            T1 = T2;
        }
    }
    for (auto &i: mu) {
        i = i/NR;
    }

    return mu;
}

std::vector<REAL> SystemWeyl::local_DOS(R3 loc) const {
    std::vector<REAL> mu = std::vector<REAL>(ORDER, 0); //weights
    //State that coraspons to that lattice site
    for (int i = 0; i < Unit_Cell_size; i++) {
        vec v(SIZE);
        v.fill(0);
        v[location(loc) + i] = 1;

        vec T0 = v;
        vec T1(SIZE);
        H(T1, v);
        vec T2(SIZE);

        std::complex<REAL> weight = 0;
        for (int j = 0; j < ORDER; j++) {
            weight = v.adjoint()*T0;
            mu[j] += weight.real()/Unit_Cell_size;
            H(T2, 2*T1);
            T2 -= T0;
            T0 = T1;
            T1 = T2;
        }
    }
    return mu;
}

std::vector<REAL> SystemWeyl::SpectralDensity(const std::function<void(vec &, const vec &)> &A) const {
    std::vector<REAL> mu = std::vector<REAL>(ORDER, 0);

    //random number engine.
    std::random_device rd;  // obtain a random seed from the OS
    std::mt19937 eng(rd());  // seed the generator
    std::uniform_real_distribution<REAL> distr(0, 2*M_PI);

#pragma omp parallel for
    for (int i = 0; i < NR; i++) {
        vec v(SIZE);
        for (auto &j: v) {
            j = std::exp(distr(eng)*COMPLEX(1i));
        }
        v.normalize();
        vec T0 = v;
        vec T0_mod = v;
        vec T1(SIZE);
        H(T1, v);
        //Weyl(T1, v);
        vec T2(SIZE);

        std::complex<REAL> weight = 0;
        for (int j = 0; j < ORDER; j++) {
            A(T0_mod, T0);
            weight = v.adjoint()*T0_mod;
            mu[j] += weight.real();
            H(T2, 2*T1);
            T2 -= T0;
            T0 = T1;
            T1 = T2;
        }
    }
    for (auto &i: mu) {
        i = i/NR;
    }

    return mu;
}

std::vector<std::vector<REAL>> SystemWeyl::SpectralDensityMulti(const h_func &A) const {
    std::vector<std::vector<REAL>> MUs;
    for (int i = 0; i < A.size(); i++) {
        MUs.emplace_back(ORDER, 0);
    }

    //random number engine.
    std::random_device rd;  // obtain a random seed from the OS
    std::mt19937 eng(rd());  // seed the generator
    std::uniform_real_distribution<REAL> distr(0, 2*M_PI);

#pragma omp parallel for
    for (int i = 0; i < NR; i++) {
        vec v(SIZE);
        for (auto &j: v) {
            j = std::exp(distr(eng)*COMPLEX(1i));
        }
        v.normalize();
        vec T0 = v;
        vec T0_mod = v;
        vec T1(SIZE);
        H(T1, v);
        vec T2(SIZE);

        std::complex<REAL> weight = 0;
        for (int j = 0; j < ORDER; j++) {
            for (int k = 0; k < MUs.size(); k++) {
                A[k](T0_mod, T0);
                weight = v.adjoint()*T0_mod;
                MUs[k][j] += weight.real();
            }

            H(T2, 2*T1);
            T2 -= T0;
            T0 = T1;
            T1 = T2;
        }
    }
    for (auto &mu: MUs) {
        for (auto &i: mu) {
            i = i/NR;
        }
    }
    return MUs;
}

std::vector<REAL>
SystemWeyl::LocalSpectralDensity(const std::function<void(vec &, const vec &)> &A, R3 loc) const {
    std::vector<REAL> mu = std::vector<REAL>(ORDER, 0); //weights
    //State that corresponds to that lattice site
    for (int i = 0; i < Unit_Cell_size; i++) {
        vec v(SIZE);
        v.fill(0);
        v[location(loc) + i] = 1;

        vec T0 = v;
        vec T0_mod = v;
        vec T1(SIZE);
        H(T1, v);
        vec T2(SIZE);

        std::complex<REAL> weight = 0;
        for (int j = 0; j < ORDER; j++) {
            A(T0_mod, T0);
            weight = v.adjoint()*T0_mod;
            mu[j] += weight.real()/Unit_Cell_size;
            H(T2, 2*T1);
            T2 -= T0;
            T0 = T1;
            T1 = T2;
        }
    }
    return mu;
}

std::vector<std::vector<REAL>> SystemWeyl::LocalSpectralDensityMulti(const h_func &A, R3 loc) const {
    std::vector<std::vector<REAL>> MUs;
    for (int i = 0; i < A.size(); i++) {
        MUs.emplace_back(ORDER, 0);
    }


#pragma omp parallel for
    for (int i = 0; i < Unit_Cell_size; i++) {
        vec v(SIZE);
        v.fill(0);
        v[location(loc) + i] = 1;

        vec T0 = v;
        vec T0_mod = v;
        vec T1(SIZE);
        H(T1, v);
        vec T2(SIZE);

        std::complex<REAL> weight = 0;
        for (int j = 0; j < ORDER; j++) {
            for (int k = 0; k < MUs.size(); k++) {
                A[k](T0_mod, T0);
                weight = v.adjoint()*T0_mod;
                MUs[k][j] += weight.real()/Unit_Cell_size;
            }

            H(T2, 2*T1);
            T2 -= T0;
            T0 = T1;
            T1 = T2;
        }
    }
    return MUs;
}

void SystemWeyl::H(vec &o, const vec &v) const {
    for(int i = 0; i < neighborsList.size(); i++) {
        Neighbor3D neighbors = neighborsList[i];
        int a = neighbors[0];
        int b = a + 1;

        COMPLEX oa = 0;
        COMPLEX ob = 0;

        //diagonal -m*sig_z term
        oa = -m*v[a];
        ob = m*v[b];

        //+x upper half Tx term
        int k = neighbors[1];
        oa += t*v[k] + ONE*tp*v[k + 1];
        ob += ONE*tp*v[k] - t*v[k + 1];
        //+x lower half
        k = neighbors[2];
        oa += t*v[k] - ONE*tp*v[k + 1];
        ob += -ONE*tp*v[k] - t*v[k + 1];


        //+y upper half Ty term
        k = neighbors[3];
        oa += t*v[k] + tp*v[k + 1];
        ob += -tp*v[k] - t*v[k + 1];
        //+y lower half
        k = neighbors[4];
        oa += t*v[k] - tp*v[k + 1];
        ob += +tp*v[k] - t*v[k + 1];


        //+z upper half
        k = neighbors[5];
        oa += t*v[k];
        ob += -t*v[k + 1];
        //+z lower half
        k = neighbors[6];
        oa += t*v[k];
        ob += -t*v[k + 1];

        //hubbard interaction
        oa += -I*M[i]*v[a];
        ob += I*M[i]*v[b];


        //scaling the energy
        o[a] = oa;
        o[b] = ob;

    }
}

/*void SystemWeyl::H(vec &o, const vec &v) const {
    for(int i = 0; i < neighborsList.size(); i++) {
        Neighbor3D neighbors = neighborsList[i];
        int a = neighbors[0];
        int b = a + 1;

        COMPLEX oa = 0;
        COMPLEX ob = 0;
        COMPLEX va = 0;
        COMPLEX vb = 0;

        //diagonal -m*sig_z term
        oa = -m*v[a];
        ob = m*v[b];

        //+x
        int plus = neighbors[1];
        int minus = neighbors[2];
        va = v[plus] + v[minus];
        vb = v[plus + 1] + v[minus + 1];
        oa += t*va + tp*vb;
        ob += tp*va - t*vb;


        //+y
        plus = neighbors[3];
        minus = neighbors[4];
        va = v[plus] + v[minus];
        vb = v[plus + 1] + v[minus + 1];
        oa += t*va - ONE*tp*vb;
        ob += ONE*tp*va - t*vb;


        //+z
        plus = neighbors[5];
        minus = neighbors[6];
        va = v[plus] + v[minus];
        vb = v[plus + 1] + v[minus + 1];
        oa += t*va;
        ob += -t*vb;


        //scaling the energy
        o[a] = oa;
        o[b] = ob;
    }
}
 */

std::vector<Neighbor3D> SystemWeyl::GenNeighbors() {
    auto list = std::vector<Neighbor3D>(Lsize);

    int counter = 0;
    for(int z = 0; z < Lz; z++) {
        for(int y = 0; y < Ly; y++){
            for(int x = 0; x < Lx; x++){

                Neighbor3D neighbor;
                neighbor[0] = location({x, y, z});

                //+x and -x
                if(x + 1 < Lx) {
                    neighbor[1] = location({x + 1, y, z});
                } else {
                    neighbor[1] = location({0, y, z});
                }

                if(x - 1 >= 0) {
                    neighbor[2] = location({x - 1, y, z});
                } else {
                    neighbor[2] = location({Lx - 1, y, z});
                }

                //+y and -y
                if(y + 1 < Ly) {
                    neighbor[3] = location({x, y + 1, z});
                } else {
                    neighbor[3] = location({x, 0, z});
                }

                if(y - 1 >= 0) {
                    neighbor[4] = location({x, y - 1, z});
                } else {
                    neighbor[4] = location({x, Ly - 1, z});
                }

                //+z and -z
                if(z + 1 < Lz) {
                    neighbor[5] = location({x, y, z + 1});
                } else {
                    neighbor[5] = location({x, y, 0});
                }

                if(z - 1 >= 0) {
                    neighbor[6] = location({x, y, z - 1});
                } else {
                    neighbor[6] = location({x, y, Lz - 1});
                }
                list[counter] = neighbor;
                counter++;
            }
        }
    }

    return list;
}

int SystemWeyl::location(R3 R) const {
    return (R[0] + R[1]*Lx + R[2]*Lx*Ly)*2;
}

int SystemWeyl::lat_location(R3 R) const {
    return (R[0] + R[1]*Lx + R[2]*Lx*Ly);
}
