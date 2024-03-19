cmake_minimum_required(VERSION 3.24)
project(KPM_4)

set(CMAKE_CXX_STANDARD 20)

find_package (Eigen3 3.4 REQUIRED)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -Wl,-ld_classic ")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wl,-ld_classic")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -Wl,-ld_classic")
    set (CPPFLAGS "/opt/homebrew/opt/libomp/include")
    set (LDFLAGS "/opt/homebrew/opt/libomp/lib")
endif()


add_executable(KPM_4 main.cpp lattice_vector.cxx lattice_vector.hxx System.cxx System.hxx Definitions.hxx SystemWeyl.cxx SystemWeyl.hxx SystemHubbard3D.cxx SystemHubbard3D.hxx SystemHubbard3D.hxx SystemHubbard2D.cxx SystemHubbard2D.hxx solver.hxx)
