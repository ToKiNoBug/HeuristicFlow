cmake_minimum_required(VERSION 3.5)

project(Heu_GA LANGUAGES CXX)

add_library(Heu_Genetic INTERFACE)
target_link_libraries(Heu_Genetic INTERFACE Heu::Global Heu::EAGlobal)

add_library(Heu::Genetic ALIAS Heu_Genetic)

file(GLOB heu_Genetic_headers "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

install(FILES ${heu_Genetic_headers} ${CMAKE_CURRENT_SOURCE_DIR}/headerCheck/InternalHeaderCheck.h
    DESTINATION include/HeuristicFlow/src/Genetic)

install(TARGETS Heu_Genetic
        DESTINATION lib
        EXPORT HeuGenetic-targets)

install(EXPORT HeuGenetic-targets
        NAMESPACE Heu::
        FILE HeuGeneticTargets.cmake
        DESTINATION lib/cmake/Heu)