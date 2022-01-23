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

#ifndef PSOBASE_HPP
#define PSOBASE_HPP

#include "PSOOption.hpp"
#include "PSOAbstrcat.hpp"

#ifndef EIGEN_CORE_H    //Detects whether libEigen is included
#ifdef OptimT_PSO_USE_EIGEN     //If user hopes to use Eigen without including it, report an error
#error You must include Eigen before you define OptimT_PSO_USE_EIGEN! Include Eigen before OptimT.
#endif
#endif

namespace OptimT {

///Abstrcat base class for most PSO solvers
template<class Var_t,size_t DIM,class Fitness_t,RecordOption Record,class...Args>
class PSOBase : public PSOAbstract<Var_t,Fitness_t,Record,Args...>
{
public:
    using Base_t = PSOAbstract<Var_t,Fitness_t,Record,Args...>;
    OPTIMT_MAKE_PSOABSTRACT_TYPES

public:
    PSOBase() {};
    virtual ~PSOBase() {};

    static size_t dimensions() {
        return DIM;
    }

protected:
    static const size_t dims=DIM;

};

///partial specialization for PSOBase with dynamic dimensions
template<class Var_t,class Fitness_t,RecordOption Record,class...Args>
class PSOBase<Var_t,Dynamic,Fitness_t,Record,Args...>
        : public PSOAbstract<Var_t,Fitness_t,Record,Args...>
{
public:
    using Base_t = PSOAbstract<Var_t,Fitness_t,Record,Args...>;
    OPTIMT_MAKE_PSOABSTRACT_TYPES

    size_t dimensions() const {
        return dims;
    }

    void setDimensions(size_t d) {
        dims=d;
    }

protected:
    size_t dims;

};

}
#endif // PSOBASE_HPP
