cmake_minimum_required(VERSION 3.5)

project(Heu_AOS LANGUAGES CXX)

add_library(Heu_AOS INTERFACE)
target_link_libraries(Heu_AOS INTERFACE Heu::Global Heu::EAGlobal)

add_library(Heu::AOS ALIAS Heu_AOS)

file(GLOB heu_AOS_headers "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

install(FILES ${heu_AOS_headers} ${CMAKE_CURRENT_SOURCE_DIR}/headerCheck/InternalHeaderCheck.h
    DESTINATION include/HeuristicFlow/src/AOS)

install(TARGETS Heu_AOS
        DESTINATION lib
        EXPORT HeuAOS-targets)

install(EXPORT HeuAOS-targets
        NAMESPACE Heu::
        FILE HeuAOSTargets.cmake
        DESTINATION lib/cmake/Heu)