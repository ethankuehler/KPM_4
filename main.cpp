#include <iostream>
#include <complex>
#include <random>
#include <iomanip>
#include <omp.h>
#include "Definitions.hxx"
#include "System.hxx"
#include "SystemWeyl.hxx"
#include "SystemHubbard3D.hxx"
#include "SystemHubbard2D.hxx"

//CONSTANT
const int NR = 12;
const int ORDER = 200;
const int Lx = 100;
const int Ly = 100;
const int Lz = 300;
const int N = ORDER;
const REAL eps = 0.1;
const REAL t = 1.0;
const REAL tp = 1.0;
const REAL m = 1.5;
const REAL I = 6;
REAL M_0 = 0.5;
REAL err = 0.01;

inline void mag(vec &o, const vec &v) {
    for (int i = 0; i < v.size(); i += 2) {
        o[i] = v[i];
        o[i + 1] = -v[i + 1];
    }
}

inline void up(vec &o, const vec &v) {
    for (int i = 0; i < v.size(); i += 2) {
        o[i] = v[i];
        o[i + 1] = 0;
    }
}

inline void down(vec &o, const vec &v) {
    for (int i = 0; i < v.size(); i += 2) {
        o[i] = 0;
        o[i + 1] = v[i + 1];
    }
}

inline void idenity(vec &o, const vec &v) {
    o = v;
}

int WriteToFile(const std::vector<REAL> &mu, const std::string &file_name) {
    // Open a CSV file in write mode
    std::ofstream outFile(file_name);

    if (!outFile) {
        // If the file couldn't be opened, print an error and exit
        std::cerr << "Failed to open the file: " << file_name << std::endl;
        return 1;
    }

    outFile << std::setprecision(std::numeric_limits<REAL>::digits10);

    // Write the array of reals to the file, each real on a new line
    for (const auto &value: mu) {
        outFile << value << ","; // Use '\n' for new line; use ',' for same line with comma separation
    }
    outFile << "0.0\n";

    // Close the file
    outFile.close();

    std::cout << "Data (" << file_name << ") written to file successfully." << std::endl;
    return 0;
}

REAL Jackson(int n, int N) {
    REAL a = (N - n + 1)*std::cos(M_PI*n/(N + 1));
    REAL b = std::sin(M_PI*n/(N + 1))*std::cos(M_PI/(N + 1))/std::sin(M_PI/(N + 1));
    return 1.0/(N + 1)*(a + b);
}

REAL Tn(REAL x, int n) {
    return std::cos(n*std::acos(x));
}

std::vector<REAL> Cheb(const std::vector<REAL> &MUs, const std::vector<REAL> &X, const std::vector<REAL> &Kernal) {
    std::vector<REAL> output(X.size());
    for (int i = 0; i < X.size(); i++) {
        REAL x = X[i];
        REAL a = 0;
        for (int n = 1; n < MUs.size(); n++) {
            REAL t = Kernal[n]*MUs[n]*Tn(x, n);
            a += t;
            if (std::isnan(t)) {
                std::cout << "fail!! kernal, mu, tn = " << Kernal[n] << ',' << MUs[n] << ',' << Tn(x, n) << '\n';
            }
        }
        output[i] = 1.0/(M_PI*std::sqrt(1 - x*x))*(Kernal[0]*MUs[0] + 2*a);
        if (std::isnan(output[i])) {
            std::cout << "meow :3\n";
            std::cout << "x = " << x << std::endl;
            std::cout << "kernal = " << Kernal[0] << std::endl;
            std::cout << "MUs = " << MUs[0] << std::endl;
            std::cout << "a = " << a << std::endl;
        }
    }
    return output;
}

REAL fermi(REAL x) {
    if (x < 0) {
        return 1;
    } else {
        return 0;
    }
}

REAL trapz(std::vector<REAL> x, std::vector<REAL> y) {
    REAL running = 0;
    for (int i = 0; i < int(x.size()) - 1; i++) {
        int a = i;
        int b = i + 1;
        running += std::abs(x[b] - x[a])/2.0*(y[b]*fermi(x[b]) + y[a]*fermi(x[a]));
    }
    return running;
}

REAL find_M_2d(std::vector<REAL> kernal, std::vector<REAL> x, REAL U, REAL err, REAL M) {
    SystemHubbard2D hubbard = SystemHubbard2D(NR, ORDER, Lx, Ly, eps, t, U, M);
    std::vector<std::function<void(vec &, const vec &)>> f;
    f.emplace_back(up);
    f.emplace_back(down);
    std::vector<std::vector<REAL>> MUs;
    std::vector<REAL> UP_DOS;
    std::vector<REAL> DOWN_DOS;
    REAL N_UP;
    REAL N_DOWN;
    REAL M_old = 0.1;
    bool stop = false;
    int i = 0;
    while (!stop) {
        std::cout << "this is the " << i << "'th iteration" << std::endl;
        MUs = hubbard.SpectralDensity(f);
        UP_DOS = Cheb(MUs[0], x, kernal);
        DOWN_DOS = Cheb(MUs[1], x, kernal);
        N_UP = trapz(x, UP_DOS);
        N_DOWN = trapz(x, DOWN_DOS);
        M = N_UP - N_DOWN;
        for (auto &k: hubbard.M) {
            k = M;
        }
        std::cout << "current U = " << U << std::endl;
        std::cout << "current M = " << M << std::endl;
        std::cout << "current <n_up> = " << N_UP << std::endl;
        std::cout << "current <n_down> = " << N_DOWN << std::endl;
        if ((std::abs(M - M_old) < err) or (std::abs(M) < err)) {
            stop = true;
        }
        M_old = M;
        i++;
    }
    return M;
}

REAL find_M_3d(std::vector<REAL> &kernal, std::vector<REAL> &x, REAL U, REAL err, REAL M) {
    SystemHubbard3D hubbard = SystemHubbard3D(NR, ORDER, Lx, Ly, Lz, eps, t, U, M);
    std::vector<REAL> mu;
    REAL M_old = M;
    bool stop = false;
    int i = 0;
    while (!stop) {
        std::cout << "this is the " << i << "'th iteration" << std::endl;
        mu = hubbard.SpectralDensity(mag);
        M = trapz(x, Cheb(mu, x, kernal));
        for (auto &k: hubbard.M) {
            k = M;
        }
        std::cout << "current U = " << U << std::endl;
        std::cout << "current M = " << M << std::endl;
        if ((std::abs(M - M_old) < err) or (std::abs(M) < err)) {
            stop = true;
        }
        M_old = M;
        i++;
    }
    return M;
}

SystemHubbard2D
find_M_2d_anti_ferro(const std::vector<REAL> &kernal, const std::vector<REAL> &x, const REAL U, const REAL err,
                     const REAL M_0) {
    SystemHubbard2D hubbard = SystemHubbard2D(NR, ORDER, Lx, Ly, eps, t, U, M_0);
    //set up anti ferromagnetic state

    //random number engine.
    std::random_device rd;  // obtain a random seed from the OS
    std::mt19937 eng(rd());  // seed the generator
    std::uniform_real_distribution<REAL> distr(-0.5, 0.5);
    REAL last = -1;

    for (auto &i: hubbard.M) {
        //i = distr(eng);
        i = last*M_0;
        last = -last;
    }

    std::vector<REAL> M(hubbard.M.size());
    std::vector<REAL> M_old(hubbard.M.size());

    M = hubbard.M;
    bool stop = false;
    int iter = 0;
    REAL AVG_OLD = 0;
    while (!stop) {
        std::cout << "this is the " << iter << "'th iteration" << std::endl;
        M_old = M;
        //calculate M at all positions
        int count = 0;
        REAL precent = 0.10;
#pragma omp parallel for shared(M_old, M, hubbard, Lx, Ly, count, std::cout, precent, kernal, x) default(none) schedule(static)
        for (int ypos = 0; ypos < Ly; ypos++) {
            for (int xpos = 0; xpos < Lx; xpos++) {
                std::vector<REAL> mu = hubbard.local_SpectralDensity(mag, {xpos, ypos});
                auto t = trapz(x, Cheb(mu, x, kernal));
                //std::cout << t << std::endl;
                M[xpos + ypos*Lx] = t;
                if (omp_get_thread_num() == 0) {
                    if (count%int(Lx*Ly/(omp_get_num_threads())*(precent)) == 0) {
                        std::cout << REAL(count)/(Lx*Ly/(omp_get_num_threads())) << " " << t << " " << std::flush;
                    }
                    count++;
                }
            }
        }
        std::cout << '\n';

        //computing root mean sqaure of M
        REAL MEAN = 0;
        for (auto &i: M) {
            MEAN += std::pow(i, 2)/M.size();
        }
        MEAN = std::sqrt(MEAN);
        std::cout << "current U = " << U << std::endl;
        std::cout << "current M RMS = " << MEAN << std::endl;

        REAL AVG = 0;
        //computing average change in M
        for (int i = 0; i < M.size(); i++) {
            AVG += std::abs(std::abs(M[i]) - std::abs(M_old[i]))/M.size();
        }
        std::cout << "current average change = " << AVG << std::endl;

        if (AVG < err or std::abs(AVG - AVG_OLD) < err) {
            std::cout << "stoping\n";
            stop = true;
        } else {
            hubbard.M = M;
            AVG_OLD = AVG;
        }
        iter++;
    }
    return hubbard;
}

std::vector<double> linspace(double start, double stop, int n) {
    std::vector<double> result;
    if (n <= 0) {
        return result; // Return an empty vector for non-positive n values
    } else if (n == 1) {
        // If n is 1, return a vector with the start value
        result.push_back(start);
        return result;
    }

    double step = (stop - start)/(n - 1);
    for (int i = 0; i < n; ++i) {
        result.push_back(start + i*step);
    }

    return result;
}

int main() {
    omp_set_num_threads(12);

    std::vector<REAL> x(N);
    for (int i = 0; i < N; i++) {
        x[i] = std::cos(M_PI*(i + 0.5)/(N));
    }
    std::reverse(x.begin(), x.end());

    x = linspace(-0.97, 0.97, N);
    if (WriteToFile(x, "../data/DOS_Plot_Points.csv") != 0) return 1;

    std::vector<REAL> kernal(ORDER);
    for (int i = 0; i < ORDER; i++) {
        kernal[i] = Jackson(i, N);
    }

    std::cout << "starting program!!!" << std::endl;

    /*
    SystemHubbard3D hubbard = SystemHubbard3D(NR, ORDER, Lx, Ly, Lz, eps, t, I, M_0);
    std::vector<REAL> M_DOS;
    REAL M;
    REAL M_old = M_0;
    */
    /*
    bool stop = false;
    int i = 0;
    while(!stop) {
        std::cout << "this is the " << i << "'th iteration" << std::endl;
        auto mu = hubbard.SpectralDensity(mag);
        auto mag_spectrum = Cheb(mu, x, kernal);
        M = trapz(x, mag_spectrum);

        std::cout << "current U = " << I << std::endl;
        std::cout << "current M = " << M << std::endl;
        if((std::abs(M-M_old) < err) or (std::abs(M) < err))  {
            stop = true;
            std::cout << "stoping\n" << '\n';
        } else {
            for(auto& k : hubbard.M) {
                k = M;
            }
        }
        M_old = M;
        i++;
    }
     */
    auto hubbard = find_M_2d_anti_ferro(kernal, x, I, err, M_0);
    //auto hubbard = SystemHubbard2D(NR, ORDER, Lx, Ly, eps, t, I, M_0);

    if (WriteToFile(hubbard.M, "../data/M.csv") != 0) return 1;

    std::cout << "done! writing output and getting final DOS\n";
    auto up_mu = hubbard.SpectralDensity(up);
    auto UP_DOS = Cheb(up_mu, x, kernal);
    if (WriteToFile(UP_DOS, "../data/UP_Plot.csv") != 0) return 1;
    auto DOWN_mu = hubbard.SpectralDensity(down);
    auto DOWN_DOS = Cheb(DOWN_mu, x, kernal);
    if (WriteToFile(DOWN_DOS, "../data/DOWN_Plot.csv") != 0) return 1;

    auto DOS = Cheb(hubbard.DOS(), x, kernal);
    if (WriteToFile(DOS, "../data/DOS_Plot.csv") != 0) return 1;
    auto N = trapz(x, DOS);
    std::cout << "Final N = " << N << std::endl;

    //DOS = Cheb(hubbard.local_DOS({Lx/2,Ly/2}), x, kernal);
    //if(WriteToFile(DOS, "../data/DOSlocal_Plot.csv") != 0) return 1;
    DOS = Cheb(hubbard.local_DOS({Lx/2 + 1, Ly/2}), x, kernal);
    if (WriteToFile(DOS, "../data/DOSlocalover_Plot.csv") != 0) return 1;

    DOS = Cheb(hubbard.local_SpectralDensity(mag, {Lx/2, Ly/2}), x, kernal);
    if (WriteToFile(DOS, "../data/DOSlocal_Plot.csv") != 0) return 1;
    N = trapz(x, DOS);
    std::cout << "Final mag  = " << N << "and M = " << hubbard.M[Lx/2 + Ly/2*Lx] << std::endl;


    DOS = Cheb(hubbard.local_SpectralDensity(mag, {Lx/2 + 1, Ly/2}), x, kernal);
    if (WriteToFile(DOS, "../data/DOSlocalover_Plot.csv") != 0) return 1;
    N = trapz(x, DOS);
    std::cout << "Final mag over = " << N << "and M = " << hubbard.M[Lx/2 + 1 + Ly/2*Lx] << std::endl;

    /*
    int N = 20;
    REAL start = 4;
    REAL end = 12;
    std::vector<REAL> Ms(N);
    std::vector<REAL> Us(N);
    for(int i = 0; i < N; i++) {
        REAL dx = (end-start)/REAL(N);
        REAL U = start;
        Ms[i] = find_M_3d(kernal, x, U + i*dx, err, M_0);
        Us[i] = U + i*dx;
    }

    if(WriteToFile(Ms, "../data/Ms_Plot.csv") != 0) return 1;
    if(WriteToFile(Us, "../data/Us_Plot.csv") != 0) return 1;
     */
}
