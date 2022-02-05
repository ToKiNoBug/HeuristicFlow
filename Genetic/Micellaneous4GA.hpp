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

#ifndef OptimT_MICELLANEOUS4GA_HPP
#define OptimT_MICELLANEOUS4GA_HPP

#include "GAAbstract.hpp"
#include "../OptimTemplates/Global"
namespace OptimT {

/*
template<typename Var_t,DivCode _min,DivCode _max>
inline void stdGAiFunNd(Var_t * v) {
    static const double min=decode<_min>::real;
    static const double max=decode<_max>::real;
    for(size_t idx=0;idx<v->size();idx++) {
        v->operator[](idx)=randD(min,max);
    }
}

template<typename Var_t,DivCode _min,DivCode _max,class Args_t>
inline void stdGAiFunNd(Var_t * v,const Args_t *) {
    stdGAiFunNd<Var_t,_min,_max>(v);
}

*/

/**
 * @brief Default crossover function for fixed-size float/double array/vector
 *        (Genetic without args)
 * 
 * @tparam Var_t type of determinate vector
 * @tparam _r crossover ratio, 0<r<1
 */
template<typename Var_t,DivCode _r>
inline void GAcFunNd(const Var_t * p1,const Var_t * p2,Var_t * c1, Var_t * c2) {
    static const double r=decode<_r>::real;
    static_assert(r>0,"r shouldn't be less than 0");
    static_assert(r<1,"r shouldn't be greater than 1");
    const size_t N=p1->size();
    for(size_t i=0;i<N;i++) {
        c1->opeartor[](i)=r*(p1->operator[](i))+(1-r)*(p2->operator[](i));
        c2->opeartor[](i)=r*(p2->operator[](i))+(1-r)*(p1->operator[](i));
    }
}

/**
 * @brief Default crossover function for dynamic-size float/double array/vector
 *        (Genetic without args)
 * 
 * @tparam Var_t type of determinate vector
 * @tparam _r crossover ratio, 0<r<1
 */
template<typename Var_t,DivCode _r>
inline void GAcFunXd(const Var_t * p1,const Var_t * p2,Var_t * c1, Var_t * c2) {
    c1->resize(p1->size());
    c2->resize(p2->size());
    GAcFunNd<Var_t,_d>(p1,p2,c1,c2);
}

/**
 * @brief Default crossover function for fixed-size float/double array/vector
 *        (Genetic with args)
 * 
 * @tparam Var_t type of determinate vector
 * @tparam _r crossover ratio, 0<r<1
 * @tparam Args_t Type of other parameter in Genetic solver
 */
template<typename Var_t,DivCode _r,class Args_t>
inline void GAcFunNd(const Var_t * p1,const Var_t * p2,Var_t * c1, Var_t * c2,const Args_t *) {
    GAcFunNd<Var_t,_r>(p1,p2,c1,c2);
}

/**
 * @brief Default crossover function for dynamic-size float/double array/vector
 *        (Genetic with args)
 * 
 * @tparam Var_t type of determinate vector
 * @tparam _r crossover ratio, 0<r<1
 * @tparam Args_t Type of other parameter in Genetic solver
 */
template<typename Var_t,DivCode _r,class Args_t>
inline void GAcFunXd(const Var_t * p1,const Var_t * p2,Var_t * c1, Var_t * c2,const Args_t *) {
    GAcFunXd<Var_t,_r>(p1,p2,c1,c2);
}


/**
 * @brief Default crossover function for fixed-size binary/symbolic array/vector
 *        (Genetic without args)
 * 
 * @tparam Var_t type of determinate vector
 */
template<typename Var_t>
inline void GAcFunNs(const Var_t * p1,const Var_t * p2,Var_t * c1, Var_t * c2) {
    const size_t N=p1->size();
    const size_t idx=randD(0,N);

    for(size_t i=0;i<N;i++) {
        if(i<idx) {
            c1->opeartor[](i)=p1->operator[](i);
            c2->opeartor[](i)=p2->operator[](i);
        }
        else {
            c1->opeartor[](i)=p2->operator[](i);
            c2->opeartor[](i)=p1->operator[](i);
        }
    }
}

/**
 * @brief Default crossover function for dynamic-size binary/symbolic array/vector
 *        (Genetic without args)
 * 
 * @tparam Var_t type of determinate vector
 */
template<typename Var_t>
inline void GAcFunXs(const Var_t * p1,const Var_t * p2,Var_t * c1, Var_t * c2) {
    c1->resize(p1->size());
    c2->resize(p2->size());
    GAcFunNs<Var_t>(p1,p2,c1,c2);
}

/**
 * @brief Default crossover function for fixed-size binary/symbolic array/vector
 *        (Genetic with args)
 * 
 * @tparam Var_t type of determinate vector
 * @tparam Args_t Type of other parameter in Genetic solver
 */
template<typename Var_t,class Args_t>
inline void GAcFunNs(const Var_t * p1,const Var_t * p2,Var_t * c1, Var_t * c2,const Args_t *) {
    GAcFunNs<Var_t>(p1,p2,c1,c2);
}

/**
 * @brief Default crossover function for dynamic-size binary/symbolic array/vector
 *        (Genetic with args)
 * 
 * @tparam Var_t type of determinate vector
 * @tparam Args_t Type of other parameter in Genetic solver
 */
template<typename Var_t,class Args_t>
inline void GAcFunXs(const Var_t * p1,const Var_t * p2,Var_t * c1, Var_t * c2,const Args_t *) {
    GAcFunXs<Var_t>(p1,p2,c1,c2);
}


}   //  namespace OptimT

#endif  //  OptimT_MICELLANEOUS4GA_HPP