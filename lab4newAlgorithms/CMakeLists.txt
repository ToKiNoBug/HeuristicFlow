cmake_minimum_required(VERSION 3.5)

#set(CMAKE_CXX_COMPILER D:/mingw/mingw64/bin/g++.exe)
#set(CMAKE_ADD2LINE D:/mingw/mingw64/bin/addr2line.exe)
#set(CMAKE_AR D:/mingw/mingw64/bin/ar.exe)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


project(OptimT_Lab LANGUAGES CXX)


include_directories(../)
include_directories(D:/CppLibs/eigen-3.4.0)

add_executable(OptimT_Lab
    lab4NSGA2.cpp
    lab4NSGA3.cpp
    makeDivCode.cpp
    main.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(OptimT_Lab PUBLIC OpenMP::OpenMP_CXX)
endif()
