cmake_minimum_required(VERSION 3.5)

project(Heu_Global LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



file(GLOB heu_Global_headers "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

install(FILES ${heu_Global_headers} ${CMAKE_CURRENT_SOURCE_DIR}/headerCheck/InternalHeaderCheck.h
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/src/Global)