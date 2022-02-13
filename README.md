# HeuristicFlow

![](https://img.shields.io/badge/C%2B%2B-17-blue?style=plastic) ![](https://img.shields.io/badge/Eigen-v3.3+-yellowgreen?style=plastic) 

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/TokiNoBug/OptimTemplates?style=plastic) ![GitHub](https://img.shields.io/github/license/TokiNoBug/OptimTemplates?style=plastic)

C++ template class lib for several optimization algorithms, like GA, NSGA2,NSGA3,PSO,etc.

HeuristicFlow doesn't depend on Eigen but it supports boosting with Eigen. Include `#include <Eigen/Dense>` **before** including OptimT will enable you to use it with Eigen.

namespace: `Heu`

Modules:
1. [Genetic](./docs/Genetic.md)
2. [Global](./docs/Genetic.md)
3. [SimpleMatrix](./docs/SimpleMatrix.md)


Implemented : 
1. [GA](./docs/Genetic/SOGA.md)
2. [NSGA-II](./docs/Genetic/NSGA2.md)
3. NSGA-III
4. PSO

Waiting to be implemented :
1. GravitySearchAlgorithm
2. BigFloodAlgorithm