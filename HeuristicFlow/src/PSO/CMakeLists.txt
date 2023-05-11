cmake_minimum_required(VERSION 3.5)

project(Heu_PSO LANGUAGES CXX)

add_library(Heu_PSO INTERFACE)
target_link_libraries(Heu_PSO INTERFACE Heu::Global Heu::EAGlobal)

add_library(Heu::PSO ALIAS Heu_PSO)

file(GLOB heu_PSO_headers "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

if(NOT ${Heu_do_install})
        return()
endif()

install(FILES ${heu_PSO_headers} ${CMAKE_CURRENT_SOURCE_DIR}/headerCheck/InternalHeaderCheck.h
        DESTINATION include/HeuristicFlow/src/PSO)

install(TARGETS Heu_PSO
        DESTINATION lib
        EXPORT HeuPSO-targets)

install(EXPORT HeuPSO-targets
        NAMESPACE Heu::
        FILE HeuPSOTargets.cmake
        DESTINATION lib/cmake/Heu)