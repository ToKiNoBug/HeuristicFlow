// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "def_TestingGenetic.h"
#include <iostream>
#include <ctime>
using namespace Heu;
using namespace std;


template<size_t M,size_t N,size_t realM=M>
void DTLZ7(const Eigen::Array<double,N,1> * x,EigenVecD_t<M> * f) {
    static_assert(realM>=2,"actural objective amount mustn't be less than 2");
    f->resize(realM,1);
    auto xm=x->template bottomRows<N-realM+1>();
    f->template topRows<realM-1>()=x->template topRows<realM-1>();
    const double g=1+9*(xm.sum())/
            //std::sqrt(1e-40+xm.square().sum())
            (N-realM+1)
            ;
    auto fi=f->template topRows<realM-1>();
    auto one_add_sum3pifi=1+(3*M_PI*fi).sin();
    const double h=realM-(fi*one_add_sum3pifi/(1+g)).sum();
    f->operator[](realM-1)=(1+g)*h;
};


void testNSGA3_DTLZ7() {
static const size_t N=20;
static const size_t M=3;
using solver_t = NSGA3<Eigen::Array<double,N,1>,
    M,
    //FITNESS_LESS_BETTER,
    DONT_RECORD_FITNESS,
    DOUBLE_LAYER,
    void>;

using Var_t = Eigen::Array<double,N,1>;
using Fitness_t = solver_t::Fitness_t;

auto iFun=GADefaults<Var_t,void,DoubleVectorOption::Eigen>::iFunNd<>;

auto DTLZ1=[](const Var_t * x,Fitness_t * f) {
    f->resize(M,1);
    auto xm=x->bottomRows<N-M+1>();
    auto xm_sub_half_square=(xm-0.5).square();
    auto cos_20pi_mul=(20*M_PI*(xm-0.5)).cos();
    const double g=100*(std::sqrt(xm.square().sum())
        +(xm_sub_half_square-cos_20pi_mul).sum());
    double accum=0.5*(1+g);
    int64_t dim=0;
    for(int64_t obj=M-1;obj>=0;obj--) {
        if(obj>0) {
            (*f)[obj]=accum*(1-(*x)[dim]);
            accum*=(*x)[dim];
            dim++;
        }
        else {
            (*f)[obj]=accum;
        }
    }
};


auto cFun=GADefaults<Var_t,void,DoubleVectorOption::Eigen>::cFunNd<encode<1,10>::code>;


auto mFun=[](const Var_t * src,Var_t * v) {
    *v=*src;
    const size_t idx=randIdx(v->size());
    double & p=v->operator[](idx);
    p+=0.5*randD(-1,1);
    p=std::min(p,1.0);
    p=std::max(p,0.0);
};

GAOption opt;
cout<<"maxGenerations=";
cin>>opt.maxGenerations;

opt.maxFailTimes=-1;
cout<<"populationSize=";
cin>>opt.populationSize;
opt.crossoverProb=0.8;
opt.mutateProb=0.1;

solver_t solver;
//solver.setObjectiveNum(M);
solver.setiFun(iFun);
solver.setfFun(DTLZ7<M,N,M>);
solver.setcFun(cFun);
solver.setmFun(mFun);
solver.setOption(opt);
solver.setReferencePointPrecision(6,4); cout<<"RPCount="<<solver.referencePointCount()<<endl;
solver.initializePop();

//cout<<"RP=["<<solver.referencePoints()<<"]';\n\n\n"<<endl;


clock_t c=clock();
solver.run();
c=clock()-c;

cout<<"solving finished in "<<c<<" ms with "<<solver.generation()<<" generations."<<endl;

cout<<"PFV=[";
for(const auto & i : solver.pfGenes()) {
    cout<<i->_Fitness.transpose()<<";\n";
}
cout<<"];\n\n\n"<<endl;

}

template<size_t M,size_t N,size_t precision>
void searchPFfun(
    size_t nIdx,
    Eigen::Array<double,N,1> * var,
    vector<pair<Eigen::Array<double,M,1>,size_t>> * dst) {

    static_assert(precision>=2,"You should assign at least 2 on a single dim");
    if(nIdx+1>=N) {
        Eigen::Array<double,M,1> f;
        for(size_t i=0;i<=precision;i++) {
            var->operator[](nIdx)=double(i)/precision;
            DTLZ7<M,N>(var,&f);
            dst->emplace_back(make_pair(f,0ULL));
        }
        return;
    }

    for(size_t i=0;i<=precision;i++) {
        var->operator[](nIdx)=double(i)/precision;
        searchPFfun<M,N,precision>(nIdx+1,var,dst);
    }

}

void searchPF() {
    static const size_t N=5,M=3;
    auto fun=DTLZ7<M,N>;
    using Var_t = Eigen::Array<double,N,1>;
    using Fitness_t = Eigen::Array<double,M,1>;
    using pair_t = pair<Fitness_t,size_t>;
    vector<pair_t> pop;
    static const size_t precision=10;
    pop.reserve(size_t(pow(N,precision)));
    Var_t v;

    cout<<"traversing"<<endl;
    clock_t c=clock();
    searchPFfun<M,N,precision>(0,&v,&pop);
    c=clock()-c;
    cout<<"traversing cost "<<c<<" ms"<<endl;

    cout<<"Nondominated sorting"<<endl;
    c=clock();
    static const int32_t thN=HfGlobal::threadNum();
#pragma omp parallel for
    for(int32_t begIdx=0;begIdx<thN;begIdx++) {
        for(int32_t idx=begIdx;idx<pop.size();idx+=thN) {
            pop[idx].second=0;
            for(size_t er=0;er<pop.size();er++) {
                pop[idx].second
                    +=internal::Pareto<M,FITNESS_LESS_BETTER>::
                    isStrongDominate(&pop[er].first,&pop[idx].first);
            }
        }
    }
    c=clock()-c;
    cout<<"NS cost "<<c<<" ms"<<endl;

    cout<<"StdPFV=[";
    for(const auto & i : pop) {
        if(i.second==0)
            cout<<i.first.transpose()<<";\n";
    }
    cout<<"];\n\n\n"<<endl;
}
