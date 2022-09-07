cmake_minimum_required(VERSION 3.5)

project(HeuristicFlow LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(EXISTS ${CMAKE_SOURCE_DIR}/compile_commands.json)
    file(REMOVE ${CMAKE_SOURCE_DIR}/compile_commands.json)
endif()


add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow)
add_subdirectory(${PROJECT_SOURCE_DIR}/test) # Test current algorithm
add_subdirectory(${PROJECT_SOURCE_DIR}/lab) # Test new algorithm


add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow/src/Genetic)
add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow/src/Global)
add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow/src/PSO)
add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow/src/EAGlobal)
add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow/src/SimpleMatrix)
add_subdirectory(${PROJECT_SOURCE_DIR}/HeuristicFlow/src/AOS)


