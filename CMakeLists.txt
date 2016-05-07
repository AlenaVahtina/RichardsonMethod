project(Richardson)
cmake_minimum_required(VERSION 2.8)
add_definitions(-std=c++11)
SET(CMAKE_BUILD_TYPE Debug)
set(SOURCE main.cpp functions.cpp plots.cpp)
set(HEADER functions.h gnuplot.h plots.h)
add_executable(${PROJECT_NAME} ${SOURCE})

