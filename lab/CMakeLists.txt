cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


project(HeuristicFlow_Lab LANGUAGES CXX)

add_definitions("-D _USE_MATH_DEFINES")

include_directories(../)
include_directories(D:/CppLibs/eigen-3.4.0)

add_executable(lab_AOS lab_AOS.cpp)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(lab_AOS PUBLIC OpenMP::OpenMP_CXX)
endif()
