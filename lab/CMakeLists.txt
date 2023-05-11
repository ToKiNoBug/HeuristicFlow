cmake_minimum_required(VERSION 3.5)

if(NOT ${Heu_build_lab})
    return()
endif()

project(HeuristicFlow_Lab LANGUAGES CXX)

add_executable(lab_AOS lab_AOS.cpp)
target_link_libraries(lab_AOS PRIVATE Heu::EAGlobal)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(lab_AOS PUBLIC OpenMP::OpenMP_CXX)
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    # add_executable(4Float128 4Float128.cpp)
    # add_compile_options(-lquadmath)
endif()