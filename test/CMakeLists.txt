cmake_minimum_required(VERSION 3.5)

if(NOT ${Heu_build_test})
    return()
endif()

project(HeuristicFlow_Test LANGUAGES CXX)

function(Heu_add_test name source dep)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE ${dep})
    add_test(NAME ${name} COMMAND ${name} --auto
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction(Heu_add_test)

Heu_add_test(Boxes Boxes.cpp Heu::Global)
Heu_add_test(TestEncoder testEncoder.cpp Heu::Global)
Heu_add_test(MinMaxCompileTime MinMaxCompileTime.cpp Heu::Global)
Heu_add_test(multiBitSet multiBitSet.cpp Heu::Global)

Heu_add_test(testFunctions testFunctions.cpp Heu::EAGlobal)

Heu_add_test(SimpleMatrix simpleMatrix.cpp Heu::SimpleMatrix)

Heu_add_test(NSGA2_Binh_And_Korn NSGA2_Binh_And_Korn.cpp Heu::Genetic)
Heu_add_test(NSGA2_Kursawe NSGA2_Kursawe.cpp Heu::Genetic)
Heu_add_test(NSGA2_ZDT3 NSGA2_ZDT3.cpp Heu::Genetic)
Heu_add_test(NSGA3_DTLZ7 NSGA3_DTLZ7.cpp Heu::Genetic)
Heu_add_test(SOGA_Ackley SOGA_Ackley.cpp Heu::Genetic)
Heu_add_test(SOGA_TSP SOGA_TSP.cpp Heu::Genetic)

Heu_add_test(PSO_RastriginFun PSO_RastriginFun.cpp Heu::PSO)
Heu_add_test(PSO_TSP PSO_TSP.cpp Heu::PSO)

Heu_add_test(AOS_Rastrigin AOS_Rastrigin.cpp Heu::AOS)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(NSGA2_Binh_And_Korn PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(NSGA2_Kursawe PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(NSGA2_ZDT3 PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(NSGA3_DTLZ7 PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(PSO_RastriginFun PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(PSO_TSP PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(SOGA_Ackley PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(SOGA_TSP PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(Boxes PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(AOS_Rastrigin PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(testFunctions PUBLIC OpenMP::OpenMP_CXX)
endif()
