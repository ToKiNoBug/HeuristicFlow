cmake_minimum_required(VERSION 3.5)

project(HeuristicFlow_Include LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(heu_include_files 
    
    ${CMAKE_CURRENT_SOURCE_DIR}/AOS 
    ${CMAKE_CURRENT_SOURCE_DIR}/EAGlobal 
    ${CMAKE_CURRENT_SOURCE_DIR}/Genetic 
    ${CMAKE_CURRENT_SOURCE_DIR}/Global 
    ${CMAKE_CURRENT_SOURCE_DIR}/PSO 
    ${CMAKE_CURRENT_SOURCE_DIR}/SimpleMatrix)

install(FILES ${heu_include_files} 
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include)