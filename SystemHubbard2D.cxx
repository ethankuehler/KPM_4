//
// Created by bruv island on 3/10/24.
//

#include "SystemHubbard2D.hxx"

SystemHubbard2D::SystemHubbard2D(int nr, int order, int lx, int ly, REAL eps, REAL t, REAL U, REAL M_0)
: NR(nr), ORDER(order), Lx(lx), Ly(ly), T(t), I(U)  {
    std::cout << "size of full Lattice: 2*" << Lx << "x" << Ly << " = " << Lx*Ly*2 << std::endl;
    Lsize = Lx * Ly;
    SIZE = Lsize*Unit_Cell_size;
    scale = (8 + I)/(2 - eps);
    Chemical_potential = 0;
    neighborsList = GenNeighbors();
    M = std::vector<REAL>(Lsize, M_0);

}

std::vector<REAL> SystemHubbard2D::DOS() const {
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

std::vector<REAL> SystemHubbard2D::local_DOS(R2 loc) const {
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

std::vector<std::vector<REAL>> SystemHubbard2D::SpectralDensity(std::vector<std::function<void(vec &, const vec &)>>& A) const {
    std::vector<std::vector<REAL>> MUs;
    for(int i = 0; i < A.size(); i++) {
        MUs.emplace_back(ORDER, 0);
    }

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
        vec T2(SIZE);

        std::complex<REAL> weight = 0;
        for (int j = 0; j < ORDER; j++) {
            for(int k = 0; k < MUs.size(); k++){
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
    for(auto& mu : MUs) {
        for (auto &i: mu) {
            i = i/NR;
        }
    }
    return MUs;
}

std::vector<REAL> SystemHubbard2D::SpectralDensity(const std::function<void(vec &, const vec &)> &A) const {
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

std::vector<REAL> SystemHubbard2D::local_SpectralDensity(const std::function<void(vec &, const vec &)> &A, R2 loc)  const{
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

void SystemHubbard2D::H(vec &o, const vec &v) const{
    for(int i = 0; i < neighborsList.size(); i++) {
        Neighbor2D neighbors = neighborsList[i];
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
        //REAL con = std::pow(M[i], 2)*(I/2);
        //o[a] += con;
        //o[b] += con;

        //chemical potential term.
        o[a] += -Chemical_potential*v[a];
        o[b] += -Chemical_potential*v[b];
        //scale energy
        o[a] = o[a]/scale;
        o[b] = o[b]/scale;
    }
}

std::vector<Neighbor2D> SystemHubbard2D::GenNeighbors() const {
    auto list = std::vector<Neighbor2D>(Lsize);
    int counter = 0;
    for(int x = 0; x < Lx; x++){
        for(int y = 0; y < Ly; y++){
            Neighbor2D neighbor;
            neighbor[0] = location({x, y});

            //+x and -x
            if(x + 1 < Lx) {
                neighbor[1] = location({x + 1, y});
            } else {
                neighbor[1] = location({0, y});
            }

            if(x - 1 >= 0) {
                neighbor[2] = location({x - 1, y});
            } else {
                neighbor[2] = location({Lx - 1, y});
            }

            //+y and -y
            if(y + 1 < Ly) {
                neighbor[3] = location({x, y + 1});
            } else {
                neighbor[3] = location({x, 0});
            }

            if(y - 1 >= 0) {
                neighbor[4] = location({x, y - 1});
            } else {
                neighbor[4] = location({x, Ly - 1});
            }

            list[counter] = neighbor;
            counter++;
        }
    }

    return list;
}

int SystemHubbard2D::location(R2 R) const {
    return (R[0] + R[1]*Ly)*2;
}
