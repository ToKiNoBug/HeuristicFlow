/*
 Copyright © 2022  TokiNoBug
This file is part of HeuristicFlow.

    HeuristicFlow is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HeuristicFlow is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HeuristicFlow.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef Heu_TYPES_HPP
#define Heu_TYPES_HPP

#include <type_traits>
#include <array>
#include <vector>
#include "Enumerations.hpp"
#include "Constants.hpp"

namespace Heu {

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


template<size_t _ObjNum>
struct initializeSize
{
    template<typename A>
    inline static void resize(A * v,size_t sz) {
        assert(v->size()==sz);
    }
};

template<>
struct initializeSize<Dynamic>
{
    template<typename A>
    inline static void resize(A * v,size_t sz) {
        v->resize(sz);
    }
};

};

#endif // Heu_TYPES_HPP
