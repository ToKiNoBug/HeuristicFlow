cmake_minimum_required(VERSION 3.5)

project(HeuristicFlow LANGUAGES CXX VERSION 1.6.3)

if(${CMAKE_SOURCE_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
    set(Heu_is_sub_project OFF)
    set(Heu_is_not_sub_project ON)
else()
    set(Heu_is_sub_project ON)
    set(Heu_is_not_sub_project OFF)
endif()

option(Heu_do_install "Install or not" ${Heu_is_not_sub_project})
option(Heu_build_test "Build tests or not" ${Heu_is_not_sub_project})
option(Heu_build_lab "Build tests or not" ${Heu_is_not_sub_project})
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Where to install")

if(NOT ${Heu_do_install})
    message(STATUS "Installation for HeuristicFlow is disabled. This function is designed for using heu as a external project.")
endif()

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

enable_testing()

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

if(${Heu_do_install})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/HeuConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/HeuConfigVersion.cmake
        DESTINATION lib/cmake/Heu)
endif()