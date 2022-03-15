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
#include "../../Global"
#include <assert.h>

namespace Heu {

namespace HeuPrivate {

template<typename Var_t,DoubleVectorOption dvo>
struct imp_GADefaults_DVO
{

    template<DivCode _r>
    inline static void imp_cFunNd(const Var_t * p1,const Var_t * p2,
                              Var_t * c1, Var_t * c2) {

        static const double constexpr r=decode<_r>::real;
        static_assert(r>0,"r shouldn't be less than 0");
        static_assert(r<1,"r shouldn't be greater than 1");

        const size_t N=p1->size();
        for(size_t i=0;i<N;i++) {
            c1->operator[](i)=r*(p1->operator[](i))+(1-r)*(p2->operator[](i));
            c2->operator[](i)=r*(p2->operator[](i))+(1-r)*(p1->operator[](i));
        }
    }

    inline static void imp_cFunSwapNs(const Var_t * p1,const Var_t * p2,
                                  Var_t * c1, Var_t * c2) {
        const size_t N=p1->size();
        const size_t idx=randD(0,N);

            c1->topRows(idx)=p1->topRows(idx);
            c2->topRows(idx)=p2->topRows(idx);
            c1->bottomRows(N-idx)=p2->bottomRows(N-idx);
            c2->bottomRows(N-idx)=p1->bottomRows(N-idx);

    }


};

template<typename Var_t>
struct imp_GADefaults_DVO<Var_t,DoubleVectorOption::Eigen>
{
    template<DivCode _r>
    inline static void imp_cFunNd(const Var_t * p1,const Var_t * p2,
                              Var_t * c1, Var_t * c2) {

        static const double constexpr r=decode<_r>::real;
        static_assert(r>0,"r shouldn't be less than 0");
        static_assert(r<1,"r shouldn't be greater than 1");

        *c1=r**p1+(1-r)**p2;
        *c2=r**p2+(1-r)**p1;
    }

    inline static void imp_cFunSwapNs(const Var_t * p1,const Var_t * p2,
                                  Var_t * c1, Var_t * c2) {
        const size_t N=p1->size();
        const size_t idx=randD(0,N);

            for(size_t i=0;i<N;i++) {
                if(i<idx) {
                    c1->operator[](i)=p1->operator[](i);
                    c2->operator[](i)=p2->operator[](i);
                }
                else {
                    c1->operator[](i)=p2->operator[](i);
                    c2->operator[](i)=p1->operator[](i);
                }
            }
    }


};

/**
 * @brief Partial specialization for GADefault struct without args
 */



}   //  namespace Heu_pri


/**
 * @brief The GADefaults struct defines several candidate operations for GA
 *
 * @tparam Var_t type of determinate vector
 * @tparam dvo double vector option of Var_t
 * @tparam Args_t type of other parameters in genetic solver
 */
template<typename Var_t,
         DoubleVectorOption dvo=DoubleVectorOption::Std,
         class Args_t=void>
struct GADefaults
{

    static_assert(!std::is_same<Args_t,void>::value,
        "The compiler run into a incorrect branch of partial specialization");

    template<BoxShape BS=SQUARE_BOX>
    inline static void iFunNd(Var_t * p,Args_t * arg) {
        static_assert(Args_t::Heu_isBox,
            "Args_t in GADefaults' paramater is not a box constraint type");
        static_assert(Args_t::encodeType==EncodeType::Real,
            "Args_t in GADefaults' paramater is real encoded");

        if constexpr (Args_t::Heu_isSquareBox) {
            for(auto & i : *p)
                i=randD(arg->min(),arg->max());
        }
        else {
            for(size_t idx=0;idx<arg->varDim();idx++)
                p->operator[](idx)=randD(arg->min()[idx],arg->max()[idx]);
        }
    }

    template<BoxShape BS=SQUARE_BOX>
    inline static void iFunNf(Var_t * p,Args_t * arg) {
        iFunNd<BS>(p,arg);
    }

    template<BoxShape BS=SQUARE_BOX>
    inline static void iFunXd(Var_t * p,Args_t * arg) {
        p->resize(arg->varDim());
        iFunNd<BS>(p,arg);
    }

    template<BoxShape BS=SQUARE_BOX>
    inline static void iFunXf(Var_t * p,Args_t * arg) {
        p->resize(arg->varDim());
        iFunNf<BS>(p,arg);
    }

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
        GADefaults<Var_t,dvo,void>::
                template cFunNd<_r>(p1,p2,c1,c2);
    }

    /**
     * @brief Default crossover function for Runtime-size float/double array/vector
     *        (Genetic with args)
     *
     * @tparam _r crossover ratio, 0<r<1
     * @tparam Args_t Type of other parameter in Genetic solver
     */
    template<DivCode _r=encode<1,5>::code>
    inline static void cFunXd(const Var_t * p1,const Var_t * p2,
                                Var_t * c1, Var_t * c2,
                                const Args_t * a) {
        GADefaults<Var_t,dvo,void>::template cFunXd<_r>(p1,p2,c1,c2);
    }

    /**
     * @brief Default crossover function for fixed-size array/vector
     *        (Genetic with args)
     *
     */
    inline static void cFunSwapNs(const Var_t * p1,const Var_t * p2,
                                Var_t * c1, Var_t * c2,
                                const Args_t *) {
        GADefaults<Var_t,dvo,void>::template cFunSwapNs(p1,p2,c1,c2);
    }

    /**
     * @brief Default crossover function for Runtime-size array/vector
     *        (Genetic with args)
     *
     */
    inline static void cFunSwapXs(const Var_t * p1,const Var_t * p2,
                                Var_t * c1, Var_t * c2,
                                const Args_t * a) {
        GADefaults<Var_t,dvo,void>::template cFunSwapXs(p1,p2,c1,c2);
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
        GADefaults<Var_t,dvo,void>::template cFunRandNs<p>(p1,p2,c1,c2);
    }

    /**
     * @brief Discrete random selection crossover by probability for Runtime-size array/vector
     *             (with args)
     *
     * @tparam p probability that c1 choose its value from p1 and c2 from p2. Default value 0.5
     */
    template<DivCode posCode=DivCode::Half>
    inline static void cFunRandXs(const Var_t * p1,const Var_t * p2,
                                  Var_t * c1, Var_t * c2,
                                  const Args_t * a) {
        GADefaults<Var_t,dvo,void>::template cFunRandXs<posCode>(p1,p2,c1,c2);
    }
};

template<typename Var_t,
         DoubleVectorOption dvo>
struct GADefaults<Var_t,dvo,void>
{

    template<DivCode _min=encode<0,1>::code,DivCode _max=encode<1,1>::code>
    inline static void iFunNd(Var_t * p) {
        static const double min=decode<_min>::real;
        static const double max=decode<_max>::real;
        //static const constexpr bool isValid=(max>min);
        //static_assert(isValid,"Max should be greater than min");

        for(size_t idx=0;idx<p->size();idx++) {
            p->operator[](idx)=randD(min,max);
        }
    }

    template<DivCode _min=encode<0,1>::code,DivCode _max=encode<1,1>::code>
    inline static void iFunNf(Var_t * p) {
        iFunNd<_min,_max>(p);
    }

    /**
     * @brief Default crossover function for fixed-size float/double array/vector
     *        (Genetic without args)
     *
     * @tparam _r crossover ratio, 0<r<1
     */
    template<DivCode _r=encode<1,5>::code>
    inline static void cFunNd(const Var_t * p1,const Var_t * p2,
                                Var_t * c1, Var_t * c2) {
        HeuPrivate::template imp_GADefaults_DVO<Var_t,dvo>::
                template imp_cFunNd<_r>(p1,p2,c1,c2);
    }

    /**
     * @brief Defaults crossover function for Runtime-size float/double array/vector
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
        HeuPrivate::template imp_GADefaults_DVO<Var_t,dvo>::
                imp_cFunSwapNs(p1,p2,c1,c2);
    }


    /**
     * @brief Default crossover function for Runtime-size array/vector
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
        static const double constexpr r=decode<p>::real;
        static_assert(r>0,"A probability shoule be greater than 0");
        static_assert(r<1,"A probability shoule be less than 1");
        const size_t N=p1->size();
        for(size_t i=0;i<N;i++) {
            c1->operator[](i)=((randD()<r)?p1:p2)->operator[](i);
            c2->operator[](i)=((randD()<r)?p2:p1)->operator[](i);
        }

    }

    /**
     * @brief Discrete random selection crossover by probability for Runtime-size array/vector
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
/*
template<typename Var_t,DoubleVectorOption dvo=DoubleVectorOption::Std,class Args_t=void>
using GADefaults =
    typename std::conditional<
    std::is_same<Args_t,void>::value,
    Heu_pri::imp_GADefaults_noParam<Var_t,dvo>,
    Heu_pri::imp_GADefaults_withParam<Var_t,dvo,Args_t>
    >::type;

*/

}   //  namespace Heu

#endif  //  Heu_MICELLANEOUS4GA_HPP
