cmake_minimum_required(VERSION 3.5)

project(HeuristicFlow_Include LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_subdirectory(./src/Genetic)
add_subdirectory(./src/Global)
add_subdirectory(./src/PSO)
add_subdirectory(./src/EAGlobal)
add_subdirectory(./src/SimpleMatrix)