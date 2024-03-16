//
// Created by bruv island on 3/6/24.
//

#include "SystemHubbard3D.hxx"

SystemHubbard3D::SystemHubbard3D(int nr, int order, int lx, int ly, int lz, REAL eps, REAL t, REAL I, REAL M_0)
: NR(nr), ORDER(order), Lx(lx), Ly(ly), Lz(lz), T(t), I(I) {
    std::cout << "size of full Lattice: 2*" << Lx << "x" << Ly << " = " << Lx*Ly*2 << std::endl;
    Lsize = Lx * Ly * Lz;
    SIZE = Lsize*Unit_Cell_size;
    scale = (12 + I)/(2 - eps);
    Chemical_potential = 0;
    neighborsList = GenNeighbors();
    M = std::vector<REAL>(Lsize, M_0);
}

std::vector<REAL> SystemHubbard3D::DOS() {
    std::vector<REAL> mu = std::vector<REAL>(ORDER, 0);

    //random number engine.
    std::random_device rd;  // obtain a random seed from the OS
    std::mt19937 eng(rd());  // seed the generator
    std::uniform_real_distribution<REAL> distr(0, 2 * M_PI);

#pragma omp parallel for
    for (int i = 0; i < NR; i++) {
        vec v(SIZE);
        for (auto &j : v){
            j = std::exp(distr(eng) * ONE);
        }
        v.normalize();
        vec T0 = v;
        vec T1(SIZE);
        H(T1, v);
        vec T2(SIZE);

        std::complex<REAL> weight = 0;
        std::stringstream ss;
        ss << "doing the " << i + 1 << " R";
        std::cout << ss.str() << std::endl;
        for (int j = 0; j < ORDER; j++) {
            weight = v.adjoint()*T0;
            mu[j] += weight.real();
            H(T2, 2*T1);
            T2 -= T0;
            T0 = T1;
            T1 = T2;
        }
    }
    for(auto& i : mu){
        i = i/NR;
    }

    return mu;
}

std::vector<REAL> SystemHubbard3D::local_DOS(R3 loc) {
    std::vector<REAL> mu = std::vector<REAL>(ORDER, 0); //weights
    //State that coraspons to that lattice site
    for(int i = 0; i < Unit_Cell_size; i++) {
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

std::vector<REAL> SystemHubbard3D::SpectralDensity(const std::function<void(vec &, const vec &)> &A) {
    std::vector<REAL> mu = std::vector<REAL>(ORDER, 0);

    //random number engine.
    std::random_device rd;  // obtain a random seed from the OS
    std::mt19937 eng(rd());  // seed the generator
    std::uniform_real_distribution<REAL> distr(0, 2 * M_PI);

#pragma omp parallel for
    for (int i = 0; i < NR; i++) {
        vec v(SIZE);
        for (auto &j : v){
            j = std::exp(distr(eng) * COMPLEX (1i));
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
    for(auto& i : mu){
        i = i/NR;
    }

    return mu;
}

std::vector<REAL> SystemHubbard3D::local_SpectralDensity(const std::function<void(vec &, const vec &)> &A, R3 loc) {
    std::vector<REAL> mu = std::vector<REAL>(ORDER, 0); //weights
    //State that coraspons to that lattice site
    for(int i = 0; i < Unit_Cell_size; i++) {
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

void SystemHubbard3D::H(vec &o, const vec &v) {
    for(int i = 0; i < neighborsList.size(); i++) {
        Neighbor3D neighbors = neighborsList[i];
        int a = neighbors[0];
        int b = a + 1;

        o[a] = 0;
        o[b] = 0;

        //tight binding
        for(int j = 1; j < neighbors.size(); j++) {
            o[a] += -T*v[neighbors[j]];
            o[b] += -T*v[neighbors[j] + 1];
        }

        //hubbard interaction
        o[a] += -I*M[i]*v[a];
        o[b] += I*M[i]*v[b];

        //constant from hubbard interaction
        //REAL con = std::pow(M[i], 2)/(2*I);
        //o[a] += con; //TODO this is needed but is off for now as it size is large
        //o[b] += con;

        //chemical potential term.
        o[a] += -Chemical_potential*v[a];
        o[b] += -Chemical_potential*v[b];
        //scale energy
        o[a] = o[a]/scale;
        o[b] = o[b]/scale;
    }
}

std::vector<Neighbor3D> SystemHubbard3D::GenNeighbors() {
    auto list = std::vector<Neighbor3D>(Lsize);

    int counter = 0;
    for(int x = 0; x < Lx; x++){
        for(int y = 0; y < Ly; y++){
            for(int z = 0; z < Lz; z++) {
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
                    neighbor[5] = location({x, y, Lz - 1});
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

int SystemHubbard3D::location(R3 R) const {
    return (R[0] + R[1]*Ly + R[2]*Lz*Ly)*2;
}
