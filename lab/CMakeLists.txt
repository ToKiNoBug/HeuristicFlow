cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


project(HeuristicFlow_Lab LANGUAGES CXX)

add_definitions("-D _USE_MATH_DEFINES")

include_directories(../)


#If you haven't added libEigen's directory into path, set the include directory manually
set(EIGEN3_INCLUDE_DIR "D:\\CppLibs\\libEigen3\\include\\eigen3")

if(NOT EIGEN3_INCLUDE_DIR)
    find_package(Eigen3 REQUIRED)
endif()

include_directories(${EIGEN3_INCLUDE_DIR})

# message("EIGEN3_INCLUDE_DIR = "${EIGEN3_INCLUDE_DIR})

add_executable(lab_AOS lab_AOS.cpp)


find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(lab_AOS PUBLIC OpenMP::OpenMP_CXX)
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    #add_executable(4Float128 4Float128.cpp)
    #add_compile_options(-lquadmath)
endif()