# HeuristicFlow

![](https://img.shields.io/badge/C%2B%2B-14-blue?style=plastic) ![](https://img.shields.io/badge/Eigen-v3.3+-yellowgreen?style=plastic) 

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/TokiNoBug/OptimTemplates?style=plastic) ![GitHub](https://img.shields.io/github/license/TokiNoBug/OptimTemplates?style=plastic)

C++ template class lib for several optimization algorithms, like GA, NSGA2,NSGA3,PSO,etc.

As a template lib, HeuristicFlow is header only.

HeuristicFlow doesn't depend on Eigen but it supports boosting with it. Include `#include <Eigen/Dense>` **before** including this lib will enable you to use it with Eigen.

HeuristicFlow has been tested on MingW(gcc 8.1.0) and MSVC143.

namespace: `Heu`

Modules:
1. [Genetic](./docs/Genetic.md)
2. [Global](./docs/Genetic.md)
3. [SimpleMatrix](./docs/SimpleMatrix.md)
4. PSO


Implemented : 
1. [GA](./docs/Genetic/SOGA.md)
2. [NSGA-II](./docs/Genetic/NSGA2.md)
3. NSGA-III
4. PSO

TODO:
1. Distinguish singular matrix better in NSGA3.
2. Simplify constructure if possible (espically for Genetic).
3. Replace `Heu::Runtime` with `Eigen::Dynamic`
4. Imporve threading implementation with OMP (learn from libEigen)
5. Change namespace.
6. Rewrite documents

Waiting to be implemented :
1. GravitySearchAlgorithm
2. BigFloodAlgorithm