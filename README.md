# HeuristicFlow

![](https://img.shields.io/badge/C%2B%2B-14-blue?style=plastic) ![](https://img.shields.io/badge/Eigen-v3.3+-yellowgreen?style=plastic) 

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/TokiNoBug/HeuristicFlow?style=plastic) ![GitHub](https://img.shields.io/github/license/TokiNoBug/HeuristicFlow?style=plastic)

C++ template class lib for several optimization algorithms, like GA, NSGA2,NSGA3,PSO,etc.

As a template lib, HeuristicFlow is header only.

Dependents: [libEigen](https://eigen.tuxfamily.org/) and **C++17** standard library.

HeuristicFlow has been tested on MinGW(gcc 8.1.0 and gcc 12.1.0), Clang(14.0.4) and MSVC(143). It works pretty well with MinGW and Clang, while may cause several warnings when with MSVC.

namespace: `heu`

Modules:
1. Genetic
2. Global
3. EAGlobal
4. SimpleMatrix
5. PSO
6. AOS


Implemented Algorithms: 
1. SOGA(Single Objective Genetic Algorithm)
2. NSGA-II(Non-dominated Sorting Genetic Algorithm II)
3. NSGA-III(Non-dominated Sorting Genetic Algorithm III)
4. PSO(Particle Swarm Optimization)
5. AOS(Atomic Orbit Search)

Implemented Testing Functions:
- Single-objective:
  - Ackley
  - Beale
  - GoldSteinPrice
  - Booth
  - Bukin
  - Matyas
  - Levy
  - Himmelblau
  - Easom
  - Cross in tray
  - Egg holder
  - Holder table
  - McCormick
  - Schaffer2
  - Schaffer4
  - Rastrigin
  - Sphere
  - Rosenbrock
  - Styblinski - Tang
- Multi-objectives:
  - Schaffer1
  - Schaffer2
  - Binh and Korn
  - Changkong Hamies
  - Poloni
  - Viennet
  - FonsecaFleming
  - Kursawe
  - DTLZ1
  - DTLZ2
  - DTLZ3
  - DTLZ4
  - DTLZ5
  - DTLZ6
  - DTLZ7
  - DTLZ8

Referenced Papers:
1. AOS<br>
   [1] [Mahdi Azizi.Atomic orbital search: A novel metaheuristic algorithm[J].Applied Mathematical Modelling.2021,93:657-693.](https://doi.org/10.1016/j.apm.2020.12.021)
2. NSGA-III<br>
   [2] [Kalyanmoy Deb,Himanshu Jain.An Evolutionary Many-Objective Optimization Algorithm Using Reference-Point-Based Nondominated Sorting Approach, Part I: Solving Problems With Box Constraint[J].IEEE Transactions on Evolutionary Computation.2014,18(4):577-601](http://dx.doi.org/10.1109/TEVC.2013.2281535)

   [3] [Himanshu Jain,Kalyanmoy Deb.An Evolutionary Many-Objective Optimization Algorithm Using Reference-point Based Non-dominated Sorting Approach, Part II: Handling Constraints and Extending to an Adaptive Approach[J].IEEE Transactions on Evolutionary Computation.2014,18(4):602-622](http://dx.doi.org/10.1109/TEVC.2013.2281534)


Waiting to be implemented :
1. GravitySearchAlgorithm
2. BigFloodAlgorithm