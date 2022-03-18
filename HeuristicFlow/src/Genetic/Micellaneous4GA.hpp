// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef Heu_MICELLANEOUS4GA_HPP
#define Heu_MICELLANEOUS4GA_HPP

#include "GAAbstract.hpp"
#include "../../Global"
#include <assert.h>

namespace Heu {

namespace internal
{

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
        const size_t idx=randIdx(N);

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

        *c1=r*(*p1)+(1-r)*(*p2);
        *c2=r*(*p2)+(1-r)*(*p1);
    }

    inline static void imp_cFunSwapNs(const Var_t * p1,const Var_t * p2,
                                  Var_t * c1, Var_t * c2) {
        const size_t N=p1->size();
        const size_t idx=randIdx(N);

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
private:
    static_assert(!std::is_same<Args_t,void>::value,
        "The compiler run into a incorrect branch of partial specialization");

    template<BoxShape BS,typename unused=void>
    struct RealBoxOp
    {
        static_assert(Args_t::Shape==BoxShape::RECTANGLE_BOX,"Wrong specialization");
        //non-square box
        inline static void imp_doiFunNd(Var_t * v,const Args_t * box) {
            for(size_t idx=0;idx<v->size();idx++) {
                v->operator[](idx)=randD(box->min()[idx],box->max()[idx]);
            }
        }

        inline static void imp_domFund_single(const Var_t * src,Var_t * v,const Args_t * box) {
            *v=*src;
            size_t idx=randIdx(v->size());
            v->operator[](idx)+=randD(-1,1)*box->learnRate()[idx];
            v->operator[](idx)=std::max(v->operator[](idx),box->min()[idx]);
            v->operator[](idx)=std::min(v->operator[](idx),box->max()[idx]);
        }

    private:
    };

    template<typename unused>
    struct RealBoxOp<BoxShape::SQUARE_BOX,unused>
    {
        //square box
        inline static void imp_doiFunNd(Var_t * v,const Args_t * box) {
            for(size_t idx=0;idx<v->size();idx++) {
                v->operator[](idx)=randD(box->min(),box->max());
            }
        }

        inline static void imp_domFund_single(const Var_t * src,Var_t * v,const Args_t * box) {
            *v=*src;
            size_t idx=randIdx(v->size());
            v->operator[](idx)+=randD(-1,1)*box->learnRate();
            v->operator[](idx)=std::max(v->operator[](idx),box->min());
            v->operator[](idx)=std::min(v->operator[](idx),box->max());
        }
    };


    template<BoxShape BS,typename unused=void>
    struct SymBoxOp //  non-square box
    {
        inline static void imp_doiFunNs(Var_t * v,const Args_t * box) {
            for(size_t idx=0;idx<v->size();idx++) {
                v->operator[](idx)=randIdx(box->min()[idx],box->max()[idx]+1);
            }
        }

        inline static void imp_domFunNs(const Var_t * src,Var_t * v,const Args_t * box) {
            *v=*src;
            const size_t idx=randIdx(v->size());
            const auto val=v->operator[](idx);
            const size_t numLess=val-box->min()[idx];
            const size_t numGreater=box->max()[idx]-val;

            if(randD()*(numLess+numGreater)<=numLess)
                v->operator[](idx)=randIdx(box->min()[idx],val);
            else
                v->operator[](idx)=randIdx(val+1,box->max()[idx]+1);
        }
    };

    template<typename unused>
    struct SymBoxOp<BoxShape::SQUARE_BOX,unused>
    {
        inline static void imp_doiFunNs(Var_t * v,const Args_t * box) {
            for(size_t idx=0;idx<v->size();idx++) {
                v->operator[](idx)=randIdx(box->min(),box->max()+1);
            }
        }

        inline static void imp_domFunNs(const Var_t * src,Var_t * v,const Args_t * box) {
            *v=*src;
            const size_t idx=randIdx(v->size());
            const auto val=v->operator[](idx);
            const size_t numLess=val-box->min();
            const size_t numGreater=box->max()-val;

            if(randD()*(numLess+numGreater)<=numLess)
                v->operator[](idx)=randIdx(box->min(),val);
            else
                v->operator[](idx)=randIdx(val+1,box->max()+1);
        }
    };

public:
    /**
     * @brief Default initialize function for fixed-sized real vectors
     *
      * @tparam unused template parameter is introduced to avoid
      * assertion failuer when non-box Args_t is used.
     */
    template<typename unused=void>
    inline static void iFunNd(Var_t * v,const Args_t * box) {
        static_assert (Args_t::isBox,
                "Default iFun requires Args_t to be a box constriant");
        static_assert (Args_t::Encoding==EncodeType::Real,
                "iFunNd requires real number encoding");
        static_assert (std::is_same<typename Args_t::Var_t,Var_t>::value,
                "Box and Var_t types must be same");

        RealBoxOp<Args_t::Shape>::imp_doiFunNd(v,box);
    }

    /**
     * @brief Default initialize function for runtime-sized real vectors
     */
    template<typename unused=void>
    inline static void iFunXd(Var_t * v,const Args_t * box) {
        v->resize(box->dimensions());
        iFunNd<unused>(v,box);
    }

    /**
     * @brief Default initialize function for fixed-sized boolean vectors
     */
    template<typename unused=void>
    inline static void iFunNb(Var_t * v,const Args_t * box) {
        static_assert (Args_t::isBox,
                "Default iFun requires Args_t to be a box constriant");
        static_assert(Args_t::Encoding==EncodeType::Binary,"iFunNb requires binary box");
        static_assert (std::is_same<typename Args_t::Var_t,Var_t>::value,
                "Box and Var_t types must be same");
        for(size_t idx=0;idx<v->size();idx++) {
            v->operator[](idx)=bool(randD()>=0.5);
        }
    }

    template<typename unused=void>
    inline static void iFunXb(Var_t * v,const Args_t * box) {
        v->resize(box->dimensions());
        iFunNb<unused>(v,box);
    }

    template<typename unused=void>
    inline static void iFunNs(Var_t * v,const Args_t * box) {
        static_assert (Args_t::isBox,
                "Default iFun requires Args_t to be a box constriant");
        static_assert(Args_t::Encoding==EncodeType::Symbolic,"iFunNs requires symbolic box");
        static_assert (std::is_same<typename Args_t::Var_t,Var_t>::value,
                "Box and Var_t types must be same");
        SymBoxOp<Args_t::Shape>::imp_doiFunNs(v,box);
    }


    template<typename unused=void>
    inline static void iFunXs(Var_t * v,const Args_t * box) {
        v->resize(box->dimensions());
        iFunNs<unused>(v,box);
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

    /**
     * @brief Default mutate function for real vectors (fixed and runtime size)
     */
    template<typename unused=void>
    inline static void mFun_d(const Var_t * src,Var_t * v,const Args_t * box) {
        static_assert (Args_t::isBox,
                "Default mFun requires Args_t to be a box constriant");
        static_assert(Args_t::Encoding==EncodeType::Real,"mFun_d requires real box");
        static_assert (std::is_same<typename Args_t::Var_t,Var_t>::value,
                "Box and Var_t types must be same");

        RealBoxOp<Args_t::Shape>::imp_domFund_single(src,v,box);
    }

    /**
     * @brief Default mutate function for binary vectors (fixed and runtime size)
     */
    template<typename unused=void>
    inline static void mFun_b(const Var_t * src,Var_t * v,const Args_t * box) {
        static_assert (Args_t::isBox,
                "Default mFun requires Args_t to be a box constriant");
        static_assert (Args_t::Encoding==EncodeType::Binary,"mFun_b requires binary box");
        static_assert (std::is_same<typename Args_t::Var_t,Var_t>::value,
                "Box and Var_t types must be same");
        *v=*src;
        size_t idx=randIdx(v->size());
        v->operator[](idx)=!v->operator[](idx);
    }

    template<typename unused=void>
    inline static void mFun_s(const Var_t * src,Var_t * v,const Args_t * box) {
        static_assert (Args_t::isBox,
                "Default mFun requires Args_t to be a box constriant");
        static_assert (Args_t::Encoding==EncodeType::Symbolic,"mFun_s requires symbolic box");
        static_assert (std::is_same<typename Args_t::Var_t,Var_t>::value,
                "Box and Var_t types must be same");

        SymBoxOp<Args_t::Shape>::imp_domFunNs(src,v,box);
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

        for(int64_t idx=0;idx<p->size();idx++) {
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
        internal::template imp_GADefaults_DVO<Var_t,dvo>::
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
        internal::template imp_GADefaults_DVO<Var_t,dvo>::
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
