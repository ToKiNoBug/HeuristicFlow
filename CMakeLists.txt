cmake_minimum_required(VERSION 3.5)

project(HeuristicFlow LANGUAGES CXX)

#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_subdirectory(./HeuristicFlow)
add_subdirectory(./test)
#add_subdirectory(./lab4newAlgorithms)
add_subdirectory(./HeuristicFlow/src/Genetic)
add_subdirectory(./HeuristicFlow/src/Global)
add_subdirectory(./HeuristicFlow/src/SimpleMatrix)
add_subdirectory(./HeuristicFlow/src/PSO)
add_subdirectory(./HeuristicFlow/src/EAGlobal)
