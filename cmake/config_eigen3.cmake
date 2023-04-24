find_package(Eigen3 3.4.0)

if(Eigen3_FOUND)
    return()
endif()

include(FetchContent)

FetchContent_Declare(Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        OVERRIDE_FIND_PACKAGE)

FetchContent_MakeAvailable(Eigen3)