cmake_minimum_required(VERSION 3.12.0)
project(lab01 LANGUAGES CXX C)

include(common/cmake-common.cmake)

add_executable(lab-09 test/lab-09.cpp test/Stokes_TI.cpp)
deal_ii_setup_target(lab-09)