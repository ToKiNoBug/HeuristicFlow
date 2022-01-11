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

#ifndef GANDS_H
#define GANDS_H
#include "./SOGA.hpp"
#include "./NSGA2.hpp"

namespace OptimT
{

template<size_t N,FitnessOption isGreaterBetter>
class GA_usual :
        public SOGA<std::array<double,N>,
                               isGreaterBetter,
                               DONT_RECORD_FITNESS,
        std::array<double,N>,std::array<double,N>,double>
{
public:
    GA_usual() {};
    ~GA_usual() {};

    using Base_t = typename SOGA<std::array<double,N>,
        FITNESS_LESS_BETTER,
        DONT_RECORD_FITNESS,
        std::array<double,N>,std::array<double,N>,double>::Base_t;
    OPTIMT_MAKE_GABASE_TYPES

    static const size_t tuple_minIdx=0;
    static const size_t tuple_maxIdx=1;
    static const size_t tuple_learningRateIdx=2;

    static void default_iFun(std::array<double,N> * x,const ArgsType * arg) {
        for(size_t idx=0;idx<N;idx++) {
            x->at(idx)=OtGlobal::randD(
                        std::get<tuple_minIdx>(*arg),
                        std::get<tuple_maxIdx>(*arg));
        }
    }

    static void default_cFun_discrete(const std::array<double,N> * p1,
                                      const std::array<double,N> * p2,
                                      std::array<double,N> * ch1,
                                      std::array<double,N> * ch2,
                                      const ArgsType *) {
        for(size_t idx=0;idx<N;idx++) {
            ch1->at(idx)=(std::rand()%2)?p1->at(idx):p2->at(idx);
            ch2->at(idx)=(std::rand()%2)?p1->at(idx):p2->at(idx);
        }
    }

    static void default_cFun_continous(const std::array<double,N> * p1,
                                      const std::array<double,N> * p2,
                                      std::array<double,N> * ch1,
                                      std::array<double,N> * ch2,
                                      const ArgsType *) {
        for(size_t idx=0;idx<N;idx++) {
            ch1->at(idx)=0.2*p1->at(idx)+(1-0.2)*p2->at(idx);
            ch2->at(idx)=0.2*p2->at(idx)+(1-0.2)*p1->at(idx);
        }
    }

    static void default_mFun(std::array<double,N> * x,const ArgsType * arg) {
        const size_t idx=size_t(OtGlobal::randD(0,N))%N;
        double & m=x->at(idx);
        m+=std::get<tuple_learningRateIdx>(*arg)*OtGlobal::randD(-1,1);
        m=std::min(m,std::get<tuple_maxIdx>(*arg)[idx]);
        m=std::max(m,std::get<tuple_minIdx>(*arg)[idx]);
    }

};

using GA2dMinimizer = GA_usual<2,FITNESS_LESS_BETTER>;
using GA3dMinimizer = GA_usual<3,FITNESS_LESS_BETTER>;
using GA4dMinimizer = GA_usual<4,FITNESS_LESS_BETTER>;
using GA5dMinimizer = GA_usual<5,FITNESS_LESS_BETTER>;
using GA6dMinimizer = GA_usual<6,FITNESS_LESS_BETTER>;

using GA2dMaximizer = GA_usual<2,FITNESS_GREATER_BETTER>;
using GA3dMaximizer = GA_usual<3,FITNESS_GREATER_BETTER>;
using GA4dMaximizer = GA_usual<4,FITNESS_GREATER_BETTER>;
using GA5dMaximizer = GA_usual<5,FITNESS_GREATER_BETTER>;
using GA6dMaximizer = GA_usual<6,FITNESS_GREATER_BETTER>;

}



#endif // GANDS_H
