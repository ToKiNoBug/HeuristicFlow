cmake_minimum_required(VERSION 3.5)
project(HeuristicFlow_Test LANGUAGES CXX)



add_executable(NSGA2_Binh_And_Korn NSGA2_Binh_And_Korn.cpp)
target_link_libraries(NSGA2_Binh_And_Korn PRIVATE Heu::Genetic)

add_executable(NSGA2_Kursawe NSGA2_Kursawe.cpp)
target_link_libraries(NSGA2_Kursawe PRIVATE Heu::Genetic)

add_executable(NSGA2_ZDT3 NSGA2_ZDT3.cpp)
target_link_libraries(NSGA2_ZDT3 PRIVATE Heu::Genetic)

add_executable(NSGA3_DTLZ7 NSGA3_DTLZ7.cpp)
target_link_libraries(NSGA3_DTLZ7 PRIVATE Heu::Genetic)

add_executable(SOGA_Ackley SOGA_Ackley.cpp)
target_link_libraries(SOGA_Ackley PRIVATE Heu::Genetic)

add_executable(SOGA_TSP SOGA_TSP.cpp)
target_link_libraries(SOGA_TSP PRIVATE Heu::Genetic)

add_executable(PSO_RastriginFun PSO_RastriginFun.cpp)
target_link_libraries(PSO_RastriginFun PRIVATE Heu::PSO)

add_executable(PSO_TSP PSO_TSP.cpp)
target_link_libraries(PSO_TSP PRIVATE Heu::PSO)

add_executable(AOS_Rastrigin AOS_Rastrigin.cpp)
target_link_libraries(AOS_Rastrigin PRIVATE Heu::AOS)

add_executable(testFunctions testFunctions.cpp)
target_link_libraries(testFunctions PRIVATE Heu::EAGlobal)

add_executable(Boxes Boxes.cpp)
target_link_libraries(Boxes PRIVATE Heu::Global)

add_executable(MinMaxCompileTime MinMaxCompileTime.cpp)
target_link_libraries(MinMaxCompileTime PRIVATE Heu::Global)

add_executable(SimpleMatrix simpleMatrix.cpp)
target_link_libraries(SimpleMatrix PRIVATE Heu::SimpleMatrix)

add_executable(TestEncoder testEncoder.cpp)
target_link_libraries(TestEncoder PRIVATE Heu::Global)

add_executable(multiBitSet multiBitSet.cpp)
target_link_libraries(multiBitSet PRIVATE Heu::Global)

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
