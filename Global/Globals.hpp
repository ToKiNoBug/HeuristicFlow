
/*
 Copyright © 2022  TokiNoBug
This file is part of Heuristic.

    Heuristic is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Heuristic is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Heuristic.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef Heu_GLOBALS_HPP
#define Heu_GLOBALS_HPP

#include <random>
#include <ctime>
#include <thread>
#include "./Enumerations.hpp"
#include "./Chaotic.hpp"
#include "./Randoms.hpp"
#include "./HeuMaths.hpp"
#include <type_traits>
#include <vector>

#include <array>

#ifndef EIGEN_CORE_H    //Detects whether libEigen is included
#ifdef Heu_GENETIC_USE_EIGEN     //If user hopes to use Eigen without including it, report an error
#error You must include Eigen before you define Heu_GENETIC_USE_EIGEN! Include Eigen before Heu.
#endif
#endif

#ifdef Heu_USE_THREADS
#include <omp.h>
#include <thread>
#endif

namespace Heu
{
///macro function for square
#define OT_square(x) (x)*(x)

///Size identifier for dynamic size (fitness or var)
const size_t Dynamic = 0;


template<typename scalar_t,size_t Dim>
using stdContainer = 
    typename std::conditional<Dim==Dynamic,
        std::vector<scalar_t>,
        std::array<scalar_t,Dim>>::type;

template<typename scalar_t,size_t Dim>
struct iniSize4StdContainer
{
public:
    inline static void iniSize(stdContainer<scalar_t,Dim> * v,size_t size) {}
};

template<typename scalar_t>
struct iniSize4StdContainer<scalar_t,Dynamic>
{

public:
    inline static void iniSize(stdContainer<scalar_t,Dynamic> * v,size_t size) {
        v->resize(size);
    }
};

#ifdef EIGEN_CORE_H
template<typename scalar_t,size_t Dim>
using EigenContainer = 
    typename std::conditional<Dim==Dynamic,
        Eigen::Array<scalar_t,Eigen::Dynamic,1>,
        Eigen::Array<scalar_t,Dim,1>>::type;  
#endif


template<size_t Size>
using stdVecD_t = stdContainer<double,Size>;

#ifdef EIGEN_CORE_H
///Array type when using Eigen array(s)
template<size_t Size>
using EigenVecD_t = EigenContainer<double,Size>;

template<DoubleVectorOption dvo,size_t Dim>
using FitnessVec_t= typename
    std::enable_if<dvo!=DoubleVectorOption::Custom,
    typename std::conditional<
    dvo==DoubleVectorOption::Eigen,
    EigenVecD_t<Dim>,
    stdVecD_t<Dim>>::type>::type;
#else
template<DoubleVectorOption dvo,size_t Dim>
using FitnessVec_t= typename
        std::enable_if<dvo==DoubleVectorOption::Std,stdVecD_t<Dim>>::type;
#endif



///Infinet value for float
const float pinfF=1.0f/0.0f;

///Infinet value for double
const double pinfD=1.0/0.0;

///negative infinet value for floa
const float nInfF=-pinfF;

///negative infinet value for double
const double ninfD=-pinfD;

///Empty class to put global variable and functions, instances of it is meanningless
class OtGlobal
{
public:
    inline static uint32_t threadNum() {
        return concurrency;
    }
private:
    OtGlobal() {};
    ~OtGlobal() {};
    static uint32_t concurrency;
};


#define Heu_MAKE_GLOBAL \
std::mt19937 Heu::global_mt19937(Heu::HeuPrivate::makeRandSeed()); \
Heu::LogisticChaos Heu::global_logistic(randD()); \
uint32_t Heu::OtGlobal::concurrency=std::thread::hardware_concurrency();

}

#endif // GLOBALS_HPP
