// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "def_TestingPSO.h"
#include <iostream>
#include <ctime>
using namespace Heu;
using namespace std;


#include <algorithm>

void testTSP_PSO(const size_t N) {
static const size_t SpaceDim=2;

    using Var_t = Eigen::ArrayXd;

    using DistanceMat_t = Eigen::ArrayXXd;

    using Solver_t = PSO_Eigen<Eigen::Dynamic,
    FITNESS_LESS_BETTER,
    RECORD_FITNESS,
    DistanceMat_t>;

    using Args_t = Solver_t::Args_t;

    Solver_t solver;

    using sortUnit = std::pair<double,size_t>;

    Solver_t::fFun_t fFun=[](const Var_t* x,const Args_t * args,double * fitness) {
        const size_t N=(*args).rows();
        std::vector<sortUnit> sortSpace(N);
        for(size_t i=0;i<N;i++) {
            sortSpace[i]=std::make_pair(x->operator[](i),i);
        }

        
        static const auto cmpFun=[](const sortUnit & a,const sortUnit & b) {
            return a.first<b.first;
        };


        std::sort(sortSpace.begin(),sortSpace.end(),cmpFun);

        *fitness=0;

        for(size_t i=0;i+1<N;i++) {
            *fitness+=(*args)(sortSpace[i].second,sortSpace[i+1].second);
        }
    };

    DistanceMat_t dMat(N,N);
    
    dMat.setZero();

    {
        Eigen::Array<double,SpaceDim,Eigen::Dynamic> points;
        points.setRandom(SpaceDim,N);

        for(size_t r=0;r<N;r++) {
            for(size_t c=0;c<N;c++) {
                if(r==c) continue;
                dMat(r,c)=(points.col(r)-points.col(c)).square().sum();
            }
        }
    }

    PSOOption opt;
    opt.inertiaFactor=0.8;
    opt.learnFactorG=2;
    opt.learnFactorP=2;
    opt.maxGeneration=50*N;
    opt.maxFailTimes=
            opt.maxGeneration/10;
    opt.populationSize=100;

    

    solver.setDimensions(N);

    solver.setPVRange(0,1,0.5);

    solver.setiFun(Solver_t::default_iFun);
    solver.setfFun(fFun);
    solver.setOption(opt);
    solver.setArgs(dMat);
    solver.initializePop();

    {
        double iniFitness;
        fFun(&solver.population().front().position,&solver.args(),&iniFitness);
        cout<<"iniFitness = "<<iniFitness<<endl;
    }

    clock_t c=clock();
    solver.run();
    c=clock()-c;

    cout<<"finished in "<<double(c)*1000/CLOCKS_PER_SEC
        <<" miliseconds and "<<solver.generation()<<" generations"<<endl;

    cout<<"result fitness = "<<solver.bestFitness()<<endl;
    /*
    cout<<"Trainning Curve=[";
    for(auto i : solver.record()) {
        cout<<i<<" , ";
    }
    cout<<"];"<<endl;
    
    cout<<"Population condition:"<<endl;
    for(const auto & i : solver.population()) {
        cout<<"fitness="<<i.fitness<<" , pBest="<<i.pBest.fitness<<endl;
    }
    */

}