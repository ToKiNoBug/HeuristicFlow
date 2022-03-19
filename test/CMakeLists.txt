cmake_minimum_required(VERSION 3.5)

#set(CMAKE_CXX_COMPILER D:/mingw/mingw64/bin/g++.exe)
#set(CMAKE_ADD2LINE D:/mingw/mingw64/bin/addr2line.exe)
#set(CMAKE_AR D:/mingw/mingw64/bin/ar.exe)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


project(HeuristicFlow_Test LANGUAGES CXX)

add_definitions("-D _USE_MATH_DEFINES")

include_directories(../)
include_directories(D:/CppLibs/eigen-3.4.0)

add_executable(HeuristicFlow_Test
    def_TestingGenetic.h
    NSGA2_Binh_And_Korn.cpp
    NSGA2_ZDT3.cpp
    SOGA_Ackley.cpp
    NSGA2_Binh_And_Korn.cpp
    NSGA3_DTLZ7.cpp
    SOGA_TSP.cpp
    def_TestingPSO.h
    PSO_RastriginFun.cpp
    PSO_testPSOBase.cpp
    PSO_TSP.cpp
    def_TestingBoxes.h
    imp_TestingBoxes.cpp
    
    main.cpp)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(HeuristicFlow_Test PUBLIC OpenMP::OpenMP_CXX)
endif()
