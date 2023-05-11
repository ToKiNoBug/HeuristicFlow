cmake_minimum_required(VERSION 3.5)

project(Heu_EAGlobal LANGUAGES CXX)

add_library(Heu_EAGlobal INTERFACE)
target_link_libraries(Heu_EAGlobal INTERFACE Heu_Global)

add_library(Heu::EAGlobal ALIAS Heu_EAGlobal)

file(GLOB heu_EAGlobal_headers "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

if(NOT ${Heu_do_install})
        return()
endif()

install(FILES ${heu_EAGlobal_headers} ${CMAKE_CURRENT_SOURCE_DIR}/headerCheck/InternalHeaderCheck.h
        DESTINATION include/HeuristicFlow/src/EAGlobal)

install(TARGETS Heu_EAGlobal
        DESTINATION lib
        EXPORT HeuEAGlobal-targets)

install(EXPORT HeuEAGlobal-targets
        NAMESPACE Heu::
        FILE HeuEAGlobalTargets.cmake
        DESTINATION lib/cmake/Heu)