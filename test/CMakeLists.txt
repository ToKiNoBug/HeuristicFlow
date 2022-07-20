cmake_minimum_required(VERSION 3.5)

#set(CMAKE_CXX_COMPILER D:/mingw/mingw64/bin/g++.exe)
#set(CMAKE_ADD2LINE D:/mingw/mingw64/bin/addr2line.exe)
#set(CMAKE_AR D:/mingw/mingw64/bin/ar.exe)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


project(HeuristicFlow_Test LANGUAGES CXX)

add_definitions("-D _USE_MATH_DEFINES")

include_directories(../)
include_directories(D:/CppLibs/eigen-3.4.0)

add_executable(NSGA2_Binh_And_Korn NSGA2_Binh_And_Korn.cpp)
add_executable(NSGA2_Kursawe NSGA2_Kursawe.cpp)
add_executable(NSGA2_ZDT3 NSGA2_ZDT3.cpp)
add_executable(NSGA3_DTLZ7 NSGA3_DTLZ7.cpp)
add_executable(PSO_RastriginFun PSO_RastriginFun.cpp)
add_executable(PSO_testPSOBase PSO_testPSOBase.cpp)
add_executable(PSO_TSP PSO_TSP.cpp)
add_executable(SOGA_Ackley SOGA_Ackley.cpp)
add_executable(SOGA_TSP SOGA_TSP.cpp)
add_executable(Boxes boxes.cpp)
add_executable(AOS_Rastrigin AOS_Rastrigin.cpp)
add_executable(testFunctions testFunctions.cpp)
add_executable(MinMaxCompileTime MinMaxCompileTime.cpp)
add_executable(SimpleMatrix simpleMatrix.cpp)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(NSGA2_Binh_And_Korn PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(NSGA2_Kursawe PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(NSGA2_ZDT3 PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(NSGA3_DTLZ7 PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(PSO_RastriginFun PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(PSO_testPSOBase PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(PSO_TSP PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(SOGA_Ackley PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(SOGA_TSP PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(Boxes PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(AOS_Rastrigin PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(testFunctions PUBLIC OpenMP::OpenMP_CXX)
endif()
