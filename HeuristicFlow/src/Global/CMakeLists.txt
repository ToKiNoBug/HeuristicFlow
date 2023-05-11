cmake_minimum_required(VERSION 3.5)

project(Heu_Global LANGUAGES CXX)

find_package(Eigen3 3.4.0 REQUIRED)

add_library(Heu_Global INTERFACE)
target_include_directories(Heu_Global INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../..>
        $<INSTALL_INTERFACE:include>)
target_link_libraries(Heu_Global INTERFACE Eigen3::Eigen)
target_compile_definitions(Heu_Global INTERFACE _USE_MATH_DEFINES)
target_compile_features(Heu_Global INTERFACE cxx_std_17)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
        target_link_libraries(Heu_Global INTERFACE OpenMP::OpenMP_CXX)
endif()

add_library(Heu::Global ALIAS Heu_Global)

file(GLOB heu_Global_headers "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

if(NOT ${Heu_do_install})
        return()
endif()

install(FILES ${heu_Global_headers} ${CMAKE_CURRENT_SOURCE_DIR}/headerCheck/InternalHeaderCheck.h
        DESTINATION include/HeuristicFlow/src/Global)

install(TARGETS Heu_Global
        DESTINATION lib
        EXPORT HeuGlobal-targets)

install(EXPORT HeuGlobal-targets
        NAMESPACE Heu::
        FILE HeuGlobalTargets.cmake
        DESTINATION lib/cmake/Heu)