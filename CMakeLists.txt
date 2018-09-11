cmake_minimum_required(VERSION 2.8)

project(ash5-cmake)
set(BOOST_PATH "C:/cpp-libs/boost_1_65_0" CACHE PATH "Path to Boost library")
#find_package(Boost REQUIRED)
set(EIGEN_PATH "C:/cpp-libs/Eigen3" CACHE PATH "Path to Eigen3 library")
#find_package(Eigen3 REQUIRED)
set(CMAKE_CXX_STANDARD 11)
include_directories("${BOOST_PATH}" "${EIGEN_PATH}")
#include_directories("${EIGEN_PATH}")
add_executable("ash5-tests" "ash5tests.cpp" "suite-config-2.ini" "suite-config-far-attractors.ini")
add_executable("ash5-suite" "ash5.cpp"
    "ash5suite.cpp" "configuration.cpp" "run_info.cpp"
    "suite-config-2.ini" "suite-config-far-attractors.ini"
    "suite-config-fix-pt-debug.ini"
    "suite-config-ns-torus.ini" "suite-config-doubled-torus.ini"
    "suite-config-chaotic.ini")
add_executable("ash5-flow" "ash5flow.cpp" "configuration.cpp" "ash5.cpp"
    "suite-config-flow.ini"
    "suite-config-ns-torus.ini" "suite-config-doubled-torus.ini"
    "suite-config-chaotic.ini")
#target_link_libraries("ash5-suite" Boost::boost)
#target_link_libraries("ash5-suite" Eigen3::Eigen)
