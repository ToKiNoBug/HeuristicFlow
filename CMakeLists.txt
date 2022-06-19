cmake_minimum_required(VERSION 3.5)

project(HeuristicFlow LANGUAGES CXX)

#set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(PROJECT_BINARY_DIR ${PROJECT_SOURCE_DIR}/../build-HeuristicFlow)
#set(CMAKE_BINARY_DIR ${PROJECT_BINARY_DIR})

#message("PROJECT_BINARY_DIR=" ${PROJECT_BINARY_DIR})
#message("CMAKE_BINARY_DIR=" ${CMAKE_BINARY_DIR})

add_subdirectory(./HeuristicFlow)
add_subdirectory(./test) # Test current algorithm
add_subdirectory(./lab) # Test new algorithm


