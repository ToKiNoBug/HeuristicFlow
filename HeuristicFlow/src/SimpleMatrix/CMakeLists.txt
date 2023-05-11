cmake_minimum_required(VERSION 3.5)

project(Heu_SimpleMatrix LANGUAGES CXX)

add_library(Heu_SimpleMatrix INTERFACE)
target_link_libraries(Heu_SimpleMatrix INTERFACE Heu::Global)

add_library(Heu::SimpleMatrix ALIAS Heu_SimpleMatrix)

file(GLOB heu_SimpleMatrix_headers "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

if(NOT ${Heu_do_install})
        return()
endif()

install(FILES ${heu_SimpleMatrix_headers} ${CMAKE_CURRENT_SOURCE_DIR}/headerCheck/InternalHeaderCheck.h
        DESTINATION include/HeuristicFlow/src/SimpleMatrix)

install(TARGETS Heu_SimpleMatrix
        DESTINATION lib
        EXPORT HeuSimpleMatrix-targets)

install(EXPORT HeuSimpleMatrix-targets
        NAMESPACE Heu::
        FILE HeuSimpleMatrixTargets.cmake
        DESTINATION lib/cmake/Heu)