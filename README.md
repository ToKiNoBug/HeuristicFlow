# HeuristicFlow

![](https://img.shields.io/badge/C%2B%2B-14-blue?style=plastic) ![](https://img.shields.io/badge/Eigen-v3.3+-yellowgreen?style=plastic) 

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/TokiNoBug/HeuristicFlow?style=plastic) ![GitHub](https://img.shields.io/github/license/TokiNoBug/HeuristicFlow?style=plastic)

C++ template class lib for several optimization algorithms, like GA, NSGA2,NSGA3,PSO,etc.

As a template lib, HeuristicFlow is header only.

HeuristicFlow doesn't depend on Eigen but it supports boosting with it. Include `#include <Eigen/Dense>` **before** including this lib will enable you to use it with Eigen.

HeuristicFlow has been tested on MingW(gcc 8.1.0) and MSVC143.

namespace: `Heu`

Modules:
1. Genetic
2. Global
3. EAGlobal
4. SimpleMatrix
5. PSO


Implemented : 
1. GA
2. NSGA-II
3. NSGA-III
4. PSO


Waiting to be implemented :
1. GravitySearchAlgorithm
2. BigFloodAlgorithm