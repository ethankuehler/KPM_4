//
// Created by bruv island on 3/2/24.
//

#include "SystemWeyl.hxx"

SystemWeyl::SystemWeyl(int nr, int order, int lx, int ly, int lz, REAL eps, REAL t, REAL tp, REAL m)
: NR(nr), ORDER(order), Lx(lx), Ly(ly), Lz(lz), t(t), tp(tp), m(m) {
    Lsize = Lx * Ly * Lz;
    SIZE = Lsize*Unit_Cell_size;
    scale = 9/(2 - eps);
    this->t = t/(2*scale);
    this->tp = tp/(2*scale);
    this->m = m/scale;
    neighborsList = GenNeighbors();
}

std::vector<REAL> SystemWeyl::DOS() {
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

std::vector<REAL> SystemWeyl::SpectralDensity(const std::function<void(vec &, const vec &)> &A) {
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
        std::cout << "doing the " << i + 1 << " R" << std::endl;
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

std::vector<Neighbor3D> SystemWeyl::GenNeighbors() {
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

int SystemWeyl::location(R3 R) const {
    return (R[0] + R[1]*Ly + R[2]*Lz*Ly)*2;
}

void SystemWeyl::H(vec &o, const vec &v) {
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
        int k = neighbors[1];;
        oa += t*v[k] + ONE*tp*v[k + 1];
        ob += ONE*tp*v[k] - t*v[k + 1];
        //+x lower half
        k = neighbors[2];;
        oa += t*v[k] - ONE*tp*v[k + 1];
        ob += -ONE*tp*v[k] - t*v[k + 1];


        //+y upper half Ty term
        k = neighbors[3];;
        oa += t*v[k] - ONE*tp*v[k + 1];
        ob += +ONE*tp*v[k] - t*v[k + 1];
        //+y lower half
        k = neighbors[4];;
        oa += t*v[k] - ONE*tp*v[k + 1];
        ob += ONE*tp*v[k] - t*v[k + 1];


        //+z upper half
        k = neighbors[5];;
        oa += t*v[k];
        ob += -t*v[k + 1];
        //+z lower half
        k = neighbors[6];;
        oa += t*v[k];
        ob += -t*v[k + 1];


        //scaling the energy
        o[a] = oa;
        o[b] = ob;

    }
}


std::vector<REAL> SystemWeyl::local_DOS(R3 loc) {
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
