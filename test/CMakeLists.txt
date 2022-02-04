cmake_minimum_required(VERSION 3.5)

#set(CMAKE_CXX_COMPILER D:/mingw/mingw64/bin/g++.exe)
#set(CMAKE_ADD2LINE D:/mingw/mingw64/bin/addr2line.exe)
#set(CMAKE_AR D:/mingw/mingw64/bin/ar.exe)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


project(OptimT_Test LANGUAGES CXX)


include_directories(../)
include_directories(D:/CppLibs/eigen-3.4.0)

add_executable(OptimT_Test
    def_TestingMatrix.h
    imp_TestingMatrix.cpp   #Test simple matrix
    def_TestingGenetic.h
    imp_TestingGenetic.cpp  #Test SOGA and NSGA2
    def_TestingPSO.h
    imp_TestingPSO.cpp  #Test PSO
    main.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(OptimT_Test PUBLIC OpenMP::OpenMP_CXX)
endif()
