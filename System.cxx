//
// Created by bruv island on 3/1/24.
//

#include "System.hxx"


TightBinding1D::TightBinding1D(REAL t, int nr, int order, int lx, REAL eps, REAL s)
: t(t), NR(nr), ORDER(order), Lx(lx){
    SIZE = Unit_Cell_size*Lx;
    this->scale = s/(2 - eps);
    gen_neighbors();

}

std::vector<REAL> TightBinding1D::DOS() {
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
        std::cout << "doing the " << i + 1 << " R" << std::endl;
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


void TightBinding1D::gen_neighbors() {
    neighbors = std::vector<Neighbor1D>(Lx);
    for(int i = 0; i < Lx; i++){
        Neighbor1D neighbor{0, 0, 0};
        neighbor[0] = i*Unit_Cell_size;
        if (i + 1 < Lx) {
            neighbor[1] = (i + 1)*Unit_Cell_size;
        } else {
            neighbor[1] = 0;
        }
        if (i - 1 > 0) {
            neighbor[2] = (i - 1)*Unit_Cell_size;
        } else {
            neighbor[2] = (Lx - 1)*Unit_Cell_size;
        }
        neighbors[i] = neighbor;
    }

}

void TightBinding1D::H(vec &out, const vec &in) {
    for(int i = 0; i < Lx; i++){
        int a = neighbors[i][0];
        int ap = neighbors[i][1];
        int am = neighbors[i][2];
        int b = a + 1;
        int bp = ap + 1;
        int bm = am + 1;
        out[a] = -t*(in[ap] + in[am])/scale;
        out[b] = -t*(in[bp] + in[bm])/scale;
    }
}

std::vector<REAL> TightBinding1D::Spectral_density(const std::function<void(vec &, const vec &)> &A) {
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
