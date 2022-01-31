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
///No-constraint N-d GA optimizer
template<size_t N,FitnessOption isGreaterBetter>
class SOGANdBox :
        public SOGA<std::array<double,N>,
                               isGreaterBetter,
                               DONT_RECORD_FITNESS,
        std::array<double,N>,std::array<double,N>,double>
{
public:
    SOGANdBox() {};
    ~SOGANdBox() {};

    using Base_t = typename SOGA<std::array<double,N>,
        FITNESS_LESS_BETTER,
        DONT_RECORD_FITNESS,
        std::array<double,N>,std::array<double,N>,double>::Base_t;
        
    OptimT_MAKE_GABASE_TYPES

    static const size_t tuple_minIdx=0;
    static const size_t tuple_maxIdx=1;
    static const size_t tuple_learningRateIdx=2;

    void setMax(const std::array<double,N> & m) {
        std::get<tuple_maxIdx>(this->args)=m;
    }

    void setMax(const double m) {
        for(auto & i : std::get<tuple_maxIdx>(this->args)) {
            i=m;
        }
    }

    void setMin(const std::array<double,N> & m) {
        std::get<tuple_minIdx>(this->args)=m;
    }
    
    void setMin(const double m) {
        for(auto & i : std::get<tuple_minIdx>(this->args)) {
            i=m;
        }
    }

    ///Default function for initializing
    static void default_iFun(std::array<double,N> * x,const ArgsType * arg) {
        for(size_t idx=0;idx<N;idx++) {
            x->at(idx)=randD(
                        std::get<tuple_minIdx>(*arg),
                        std::get<tuple_maxIdx>(*arg));
        }
    }
    ///Default function for crossover
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
    ///Default function for crossover
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
    ///Default function for mutation
    static void default_mFun(std::array<double,N> * x,const ArgsType * arg) {
        const size_t idx=size_t(randD(0,N))%N;
        double & m=x->at(idx);
        m+=std::get<tuple_learningRateIdx>(*arg)*randD(-1,1);
        m=std::min(m,std::get<tuple_maxIdx>(*arg)[idx]);
        m=std::max(m,std::get<tuple_minIdx>(*arg)[idx]);
    }
private:
    static_assert(N>=1,"Dimensional numbers less than 1!");

};


}



#endif // GANDS_H
