
/*
 Copyright © 2022  TokiNoBug
This file is part of OptimTemplates.

    OptimTemplates is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OptimTemplates is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OptimTemplates.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef OptimT_GLOBALS_HPP
#define OptimT_GLOBALS_HPP

#include <random>
#include <ctime>
#include <thread>
#include "./Enumerations.hpp"
#include "./Chaotic.hpp"
#include "./Randoms.hpp"
#include "./OptimTMaths.hpp"
#include <type_traits>
#include <vector>

#include <array>

#ifdef OptimT_DO_PARALLELIZE
#include <omp.h>
#include <thread>
#endif

namespace OptimT
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
    inline static void iniSize(stdContainer<scalar_t,Dim> * v,size_t size) {
        v->resize(size);
    }
};

#ifdef EIGEN_CORE_H
template<typename scalar_t,size_t Dim>
using EigenContainer = 
    typename std::conditional<Dim==Dynamic,
        Eigen::ArrayXd<scalar_t>,
        Eigen::Array<scalar_t,Dim,1>>::type;  
#endif


template<size_t Size>
using stdVecD_t = typename stdContainer<double,Size>;

#ifdef EIGEN_CORE_H
///Array type when using Eigen array(s)
template<size_t Size>
using EigenVecD_t = typename EigenContainer<double,Size>;

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


#define OptimT_MAKE_GLOBAL \
std::mt19937 OptimT::global_mt19937(OptimT::OptimTPrivate::makeRandSeed()); \
OptimT::LogisticChaos OptimT::global_logistic(randD()); \
uint32_t OptimT::OtGlobal::concurrency=std::thread::hardware_concurrency();

}

#endif // GLOBALS_HPP
