cmake_minimum_required(VERSION 3.5)

project(OptimT_Test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(../)
include_directories(D:/CppLibs/eigen-3.4.0)

add_executable(OptimT_Test
    testNsga2.cpp   #Manual implementation of NSGA2 inheriting GABase
    imp_TestingMatrix.cpp   #Test simple matrix
    imp_TestingGenetic.cpp  #Test SOGA and NSGA2
    imp_TestingPSO.cpp  #Test PSO
    main.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(OptimT_Test PUBLIC OpenMP::OpenMP_CXX)
endif()
