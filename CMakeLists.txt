cmake_minimum_required(VERSION 3.5)

project(HeuristicFlow LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(PROJECT_BINARY_DIR ${PROJECT_SOURCE_DIR}/.${PROJECT_SOURCE_DIR}/build-HeuristicFlow)
#set(CMAKE_BINARY_DIR ${PROJECT_BINARY_DIR})

#message("PROJECT_BINARY_DIR=" ${PROJECT_BINARY_DIR})
#message("CMAKE_BINARY_DIR=" ${CMAKE_BINARY_DIR})

add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow)
add_subdirectory(${PROJECT_SOURCE_DIR}/test) # Test current algorithm
add_subdirectory(${PROJECT_SOURCE_DIR}/lab) # Test new algorithm


add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow/src/Genetic)
add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow/src/Global)
add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow/src/PSO)
add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow/src/EAGlobal)
add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow/src/SimpleMatrix)
add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow/src/AOS)

#set(PROJECT_BINARY_DIR ${PROJECT_SOURCE_DIR}/../build-HeuristicFlow)
#set(CMAKE_BINARY_DIR ${PROJECT_BINARY_DIR})

#message("PROJECT_BINARY_DIR=" ${PROJECT_BINARY_DIR})
#message("CMAKE_BINARY_DIR=" ${CMAKE_BINARY_DIR})

add_subdirectory(./HeuristicFlow)
add_subdirectory(./test)
add_subdirectory(./lab)
#add_subdirectory(./lab4newAlgorithms)
add_subdirectory(./HeuristicFlow/src/Genetic)
add_subdirectory(./HeuristicFlow/src/Global)
add_subdirectory(./HeuristicFlow/src/PSO)
add_subdirectory(./HeuristicFlow/src/EAGlobal)

