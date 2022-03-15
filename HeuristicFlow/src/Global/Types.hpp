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
    typename std::conditional<Dim==Runtime,
        std::vector<scalar_t>,
        std::array<scalar_t,Dim>>::type;

#ifdef EIGEN_CORE_H
template<typename scalar_t,size_t Dim>
using EigenContainer =
    typename std::conditional<Dim==Runtime,
        Eigen::Array<scalar_t,Eigen::Dynamic,1>,
        Eigen::Array<scalar_t,Dim,1>>::type;
#endif


template<size_t Size>
using stdVecD_t = stdContainer<double,Size>;


#ifdef EIGEN_CORE_H
template<typename scalar_t,size_t Size,DoubleVectorOption DVO>
using Container = typename std::conditional<
    (DVO!=DoubleVectorOption::Eigen),
    stdContainer<scalar_t,Size>,
    EigenContainer<scalar_t,Size>>::type;
#else
template<typename scalar_t,size_t Size,DoubleVectorOption DVO>
using Container = typename std::enable_if
    <(DVO!=DoubleVectorOption::Eigen),
    stdContainer<scalar_t,Size>>::type;
#endif

#ifdef EIGEN_CORE_H
///Array type when using Eigen array(s)
template<size_t Size>
using EigenVecD_t = EigenContainer<double,Size>;
#endif

template<DoubleVectorOption dvo,size_t Dim>
using FitnessVec_t= Container<double,Dim,dvo>;

template<size_t _ObjNum>
struct initializeSize
{
    template<typename A>
    inline static void resize(A * v,size_t sz) {
        assert(v->size()==sz);
    }
};

template<>
struct initializeSize<Runtime>
{
    template<typename A>
    inline static void resize(A * v,size_t sz) {
        v->resize(sz);
    }
};

};

#endif // Heu_TYPES_HPP
