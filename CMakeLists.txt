cmake_minimum_required(VERSION 3.5)

project(HeuristicFlow LANGUAGES CXX VERSION 1.6.2.1)

set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Where to install")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(EXISTS ${CMAKE_SOURCE_DIR}/compile_commands.json)
    file(REMOVE ${CMAKE_SOURCE_DIR}/compile_commands.json)
endif()

set(Heu_submodules
        Global
        SimpleMatrix
        EAGlobal
        AOS
        Genetic
        PSO)


if(NOT DEFINED EIGEN3_INCLUDE_DIR)
    include(cmake/config_eigen3.cmake)
endif()
message(STATUS "EIGEN3_INCLUDE_DIR = " ${EIGEN3_INCLUDE_DIR})


add_subdirectory(HeuristicFlow)
add_subdirectory(test) # Test current algorithm
add_subdirectory(lab) # Test new algorithm

include(CMakePackageConfigHelpers)
configure_package_config_file(HeuConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/HeuConfig.cmake
        INSTALL_DESTINATION lib/cmake/Heu)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/HeuConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/HeuConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/HeuConfigVersion.cmake
        DESTINATION lib/cmake/Heu)