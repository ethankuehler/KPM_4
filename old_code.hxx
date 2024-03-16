#include <iostream>
#include <complex>
#include <random>
#include <iomanip>
#include "Definitions.hxx"
#include "System.hxx"
#include "SystemWeyl.hxx"
#include "SystemHubbard3D.hxx"
#include "SystemHubbard2D.hxx"
#include "main.cpp"

#ifndef KPM_4_OLD_CODE_HXX
#define KPM_4_OLD_CODE_HXX

int old(){
    //SystemWeyl weyl = SystemWeyl(NR, ORDER, Lx, Lx, Lz, eps, t, tp, m);
    //SystemHubbard2D hubbard = SystemHubbard2D(NR, ORDER, Lx, Ly, eps, t, 8.1);
    std::cout << "starting program!!!" << std::endl;
    std::cout << "size of full Lattice: 2*" << Lx << "x" << Ly << "*" << Lz << " = " << Lx*Ly*2 << std::endl;
    std::cout << "finding MU_DOS" << std::endl;

    /*std::vector<std::function<void(vec &, const vec &)>> f;
    f.emplace_back(idenity);
    f.emplace_back(up);
    f.emplace_back(down);
    auto MUs = hubbard.SpectralDensity(f);
    //std::cout << "finding MU_MAG" << std::endl;
    //std::vector<REAL> mu_mag = weyl.SpectralDensity(mag);
    if(WriteToFile(MUs[0], "../data/mu_dos.csv") != 0) return 1;
    //if(WriteToFile(mu_mag, "../data/mu_mag.csv") != 0) return 1;

    std::vector<REAL> x(2*ORDER);
    for(int i = 0; i < 2*ORDER; i++){
        x[i] = std::cos(M_PI*(i + 0.5)/(2*ORDER));
    }
    if(WriteToFile(x, "../data/DOS_Plot_Points.csv") != 0) return 1;

    std::vector<REAL> kernal(ORDER);
    for(int i = 0; i < ORDER; i++) {
        kernal[i] = Jackson(i, ORDER);
    }

    auto DOS = Cheb(MUs[0], x, kernal);
    if(WriteToFile(DOS, "../data/DOS_Plot.csv") != 0) return 1;

    auto UP_DOS = Cheb(MUs[1], x, kernal);
    if(WriteToFile(UP_DOS, "../data/UP_Plot.csv") != 0) return 1;
    auto DOWN_DOS = Cheb(MUs[2], x, kernal);
    if(WriteToFile(DOWN_DOS, "../data/DOWN_Plot.csv") != 0) return 1;
    */

    //ok lets try to self-consit

    std::vector<REAL> x(2*ORDER);
    for(int i = 0; i < 2*ORDER; i++){
        x[i] = std::cos(M_PI*(i + 0.5)/(2*ORDER));
    }
    if(WriteToFile(x, "../data/DOS_Plot_Points.csv") != 0) return 1;

    std::vector<REAL> kernal(ORDER);
    for(int i = 0; i < ORDER; i++) {
        kernal[i] = Jackson(i, ORDER);
    }


    SystemHubbard2D hubbard = SystemHubbard2D(NR, ORDER, Lx, Ly, eps, t, U, M_0);
    std::vector<std::function<void(vec &, const vec &)>> f;
    f.emplace_back(up);
    f.emplace_back(down);
    f.emplace_back(mag);
    std::vector<std::vector<REAL>> MUs;
    std::vector<REAL> UP_DOS;
    std::vector<REAL> DOWN_DOS;
    std::vector<REAL> mag;
    REAL N_UP;
    REAL N_DOWN;
    REAL M;
    REAL M_old = M_0;

    bool stop = false;
    int i = 0;
    while(!stop) {
        std::cout << "this is the " << i << "'th iteration" << std::endl;
        MUs = hubbard.SpectralDensity(f);
        UP_DOS = Cheb(MUs[0], x, kernal);
        DOWN_DOS = Cheb(MUs[1], x, kernal);
        mag = Cheb(MUs[0], x, kernal);
        N_UP = trapz(x, UP_DOS);
        N_DOWN = trapz(x, DOWN_DOS);
        auto m = trapz(x, mag);
        M = N_UP - N_DOWN;
        for(auto& k : hubbard.M) {
            k = M;
        }
        std::cout << "current U = " << U << std::endl;
        std::cout << "current M = " << M << std::endl;
        std::cout << "M compared to m = " << M << " " << m << '\n';
        std::cout << "current <n_up> = " << N_UP << std::endl;
        std::cout << "current <n_down> = " << N_DOWN << std::endl;
        if((std::abs(M-M_old) < err) or (std::abs(M) < err))  {
            stop = true;
        }
        M_old = M;
        i++;
    }
    std::cout << "done! writing output and getting final DOS";
    if(WriteToFile(UP_DOS, "../data/UP_Plot.csv") != 0) return 1;
    if(WriteToFile(DOWN_DOS, "../data/DOWN_Plot.csv") != 0) return 1;

    auto DOS = Cheb(hubbard.DOS(), x, kernal);
    auto N = trapz(x, DOS);
    std::cout << "Final N = " << N << std::endl;

    if(WriteToFile(DOS, "../data/DOS_Plot.csv") != 0) return 1;

    /*
    int N = 20;
    REAL start = 2.5;
    REAL end = 3.5;
    std::vector<REAL> Ms(N);
    std::vector<REAL> Us(N);
    for(int i = 0; i < N; i++) {
        REAL dx = (end-start)/REAL(N);
        REAL U = start;
        Ms[i] = find_M(kernal, x, U + i*dx, 0.001);
        Us[i] = U + i*dx;
    }

    if(WriteToFile(Ms, "../data/Ms_Plot.csv") != 0) return 1;
    if(WriteToFile(Us, "../data/Us_Plot.csv") != 0) return 1;

    */

    return 0;
}

#endif //KPM_4_OLD_CODE_HXX
