project(Richardson)
cmake_minimum_required(VERSION 2.8)
add_definitions(-std=c++11)
SET(CMAKE_BUILD_TYPE Debug)
set(SOURCE main.cpp plots.cpp matrix/basematrix.cpp matrix/crsmatrix.cpp matrix/normalmatrix.cpp
    richardsonMethodWithChebyshevOrderedSetOfParameters.cpp common.cpp)
set(HEADER gnuplot.h plots.h matrix/basematrix.h matrix/crsmatrix.h matrix/normalmatrix.h
    richardsonMethodWithChebyshevOrderedSetOfParameters.h common.h)
add_executable(${PROJECT_NAME} ${SOURCE} ${HADER})

