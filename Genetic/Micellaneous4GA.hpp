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

#ifndef Heu_MICELLANEOUS4GA_HPP
#define Heu_MICELLANEOUS4GA_HPP

#include "GAAbstract.hpp"
#include "../HeuristicFlow/Global"
#include <assert.h>

namespace Heu {

/*
template<DivCode _min,DivCode _max>
inline void stdGAiFunNd(Var_t * v) {
    static const double min=decode<_min>::real;
    static const double max=decode<_max>::real;
    for(size_t idx=0;idx<v->size();idx++) {
        v->operator[](idx)=randD(min,max);
    }
}

template<DivCode _min,DivCode _max>
inline void stdGAiFunNd(Var_t * v,const Args_t *) {
    stdGAiFunNd<_min,_max>(v);
}

*/

namespace Heu_pri
{

/**
 * @brief Partial specialization for GADefault struct without args
 */
template<typename Var_t,DoubleVectorOption dvo>
struct imp_GADefaults_noParam
{
    /**
     * @brief Default crossover function for fixed-size float/double array/vector
     *        (Genetic without args)
     *
     * @tparam _r crossover ratio, 0<r<1
     */
    template<DivCode _r=encode<1,5>::code>
    inline static void cFunNd(const Var_t * p1,const Var_t * p2,
                                Var_t * c1, Var_t * c2) {
#define Heu_PRIVATE_IMP_cFunNd \
        static const double constexpr r=decode<_r>::real; \
        static_assert(r>0,"r shouldn't be less than 0"); \
        static_assert(r<1,"r shouldn't be greater than 1"); \
        if constexpr (dvo==DoubleVectorOption::Eigen) { \
            *c1=r**p1+(1-r)**p2; \
            *c2=r**p2+(1-r)**p1; \
        } \
        else { \
            const size_t N=p1->size(); \
            for(size_t i=0;i<N;i++) { \
            c1->operator[](i)=r*(p1->operator[](i))+(1-r)*(p2->operator[](i)); \
            c2->operator[](i)=r*(p2->operator[](i))+(1-r)*(p1->operator[](i)); \
        } \
        }
        

        Heu_PRIVATE_IMP_cFunNd
    }

    /**
     * @brief Defaults crossover function for dynamic-size float/double array/vector
     *        (Genetic without args)
     *
     * @tparam _r crossover ratio, 0<r<1
     */
    template<DivCode _r=encode<1,5>::code>
    inline static void cFunXd(const Var_t * p1,const Var_t * p2,
                                Var_t * c1, Var_t * c2) {
#define Heu_PRIVATE_IMP_cFunX \
        c1->resize(p1->size()); \
        c2->resize(p2->size());

        Heu_PRIVATE_IMP_cFunX

        cFunNd<_r>(p1,p2,c1,c2);
    }

    /**
     * @brief Default crossover function for fixed-size array/vector
     *        (Genetic without args)
     *
     */

    inline static void cFunSwapNs(const Var_t * p1,const Var_t * p2,
                                Var_t * c1, Var_t * c2) {
#define Heu_PRIVATE_IMP_cFunSwapNs \
        const size_t N=p1->size(); \
        const size_t idx=randD(0,N); \
        if constexpr (dvo==DoubleVectorOption::Eigen) { \
            c1->topRows(idx)=p1->topRows(idx); \
            c2->topRows(idx)=p2->topRows(idx); \
            c1->bottomRows(N-idx)=p2->bottomRows(N-idx); \
            c2->bottomRows(N-idx)=p1->bottomRows(N-idx); \
        } \
        else { \
            for(size_t i=0;i<N;i++) { \
                if(i<idx) { \
                    c1->operator[](i)=p1->operator[](i); \
                    c2->operator[](i)=p2->operator[](i); \
                } \
                else { \
                    c1->operator[](i)=p2->operator[](i); \
                    c2->operator[](i)=p1->operator[](i); \
                } \
            } \
        }

        Heu_PRIVATE_IMP_cFunSwapNs

    }


    /**
     * @brief Default crossover function for dynamic-size array/vector
     *        (Genetic without args)
     *
     */
    inline static void cFunSwapXs(const Var_t * p1,const Var_t * p2,
                                Var_t * c1, Var_t * c2) {

        Heu_PRIVATE_IMP_cFunX

        cFunSwapNs(p1,p2,c1,c2);
    }

    /**
     * @brief Discrete random selection crossover by probability for fixed-size array/vector
     *             (without args)
     *
     * @tparam p probability that c1 choose its value from p1 and c2 from p2. Default value 0.5
     */
    template<DivCode p=DivCode::Half>
    inline static void cFunRandNs(const Var_t * p1,const Var_t * p2,
                                  Var_t * c1, Var_t * c2) {
#define Heu_PRIVATE_IMP_cFunRandNs \
        static const double constexpr r=decode<p>::real; \
        static_assert(r>0,"A probability shoule be greater than 0"); \
        static_assert(r<1,"A probability shoule be less than 1"); \
        const size_t N=p1->size(); \
        for(size_t i=0;i<N;i++) { \
            c1->operator[](i)=((randD()<r)?p1:p2)->operator[](i); \
            c2->operator[](i)=((randD()<r)?p2:p1)->operator[](i); \
        }

        Heu_PRIVATE_IMP_cFunRandNs

    }

    /**
     * @brief Discrete random selection crossover by probability for dynamic-size array/vector
     *             (without args)
     *
     * @tparam p probability that c1 choose its value from p1 and c2 from p2. Default value 0.5
     */
    template<DivCode p=DivCode::Half>
    inline static void cFunRandXs(const Var_t * p1,const Var_t * p2,
                                  Var_t * c1, Var_t * c2) {

        Heu_PRIVATE_IMP_cFunX

        cFunRandNs(p1,p2,c1,c2);
    }
};


/**
 * @brief The GADefaults struct defines several candidate operations for GA
 *
 * @tparam Var_t type of determinate vector
 * @tparam Args_t type of other parameters in genetic solver
 */

template<typename Var_t,class Args_t,DoubleVectorOption dvo>
struct imp_GADefaults_withParam
{
    static_assert(!std::is_same<Args_t,void>::value,
        "The compiler run into a incorrect branch of partial specialization");



    /**
     * @brief Default crossover function for fixed-size float/double array/vector
     *        (Genetic with args)
     *
     * @tparam _r crossover ratio, 0<r<1
     * @tparam Args_t Type of other parameter in Genetic solver
     */
    template<DivCode _r=encode<1,5>::code>
    inline static void cFunNd(const Var_t * p1,const Var_t * p2,
                                Var_t * c1, Var_t * c2,
                                const Args_t *) {
        Heu_PRIVATE_IMP_cFunNd
    }

    /**
     * @brief Default crossover function for dynamic-size float/double array/vector
     *        (Genetic with args)
     *
     * @tparam _r crossover ratio, 0<r<1
     * @tparam Args_t Type of other parameter in Genetic solver
     */
    template<DivCode _r=encode<1,5>::code>
    inline static void cFunXd(const Var_t * p1,const Var_t * p2,
                                Var_t * c1, Var_t * c2,
                                const Args_t * a) {
        Heu_PRIVATE_IMP_cFunX
        cFunNd<_r>(p1,p2,c1,c2,a);
    }

    /**
     * @brief Default crossover function for fixed-size array/vector
     *        (Genetic with args)
     *
     */
    inline static void cFunSwapNs(const Var_t * p1,const Var_t * p2,
                                Var_t * c1, Var_t * c2,
                                const Args_t *) {
        Heu_PRIVATE_IMP_cFunSwapNs
    }

    /**
     * @brief Default crossover function for dynamic-size array/vector
     *        (Genetic with args)
     *
     */
    inline static void cFunSwapXs(const Var_t * p1,const Var_t * p2,
                                Var_t * c1, Var_t * c2,
                                const Args_t * a) {
        Heu_PRIVATE_IMP_cFunX
        cFunSwapXs(p1,p2,c1,c2,a);
    }


    /**
     * @brief Discrete random selection crossover by probability for fixed-size array/vector
     *             (with args)
     *
     * @tparam p probability that c1 choose its value from p1 and c2 from p2. Default value 0.5
     */
    template<DivCode p=DivCode::Half>
    inline static void cFunRandNs(const Var_t * p1,const Var_t * p2,
                                  Var_t * c1, Var_t * c2,
                                  const Args_t *) {
        Heu_PRIVATE_IMP_cFunRandNs
    }

    /**
     * @brief Discrete random selection crossover by probability for dynamic-size array/vector
     *             (with args)
     *
     * @tparam p probability that c1 choose its value from p1 and c2 from p2. Default value 0.5
     */
    template<DivCode posCode=DivCode::Half>
    inline static void cFunRandXs(const Var_t * p1,const Var_t * p2,
                                  Var_t * c1, Var_t * c2,
                                  const Args_t * a) {
        Heu_PRIVATE_IMP_cFunX

        cFunRandNs<posCode>(p1,p2,c1,c2,a);
    }

};

}   //  namespace Heu_pri



template<typename Var_t,DoubleVectorOption dvo=DoubleVectorOption::Std,class Args_t=void>
using GADefaults =
    typename std::conditional<
    std::is_same<Args_t,void>::value,
    Heu_pri::imp_GADefaults_noParam<Var_t,dvo>,
    Heu_pri::imp_GADefaults_withParam<Var_t,Args_t,dvo>
    >::type;


}   //  namespace Heu

#endif  //  Heu_MICELLANEOUS4GA_HPP
