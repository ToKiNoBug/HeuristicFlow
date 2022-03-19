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


void testRastriginFun() {
    static const size_t N=20;
    using solver_t = PSO_Eigen<N,
    FITNESS_LESS_BETTER,
    RECORD_FITNESS>;
    
    using Var_t = EigenVecD_t<N>;

    PSOOption opt;
    opt.populationSize=400;
    opt.maxGeneration=50*N;
    opt.maxFailTimes=
            -1;
            //opt.maxGeneration/10;
    opt.inertiaFactor=0.8;
    opt.learnFactorG=2;
    opt.learnFactorP=2;


    solver_t solver;

    solver_t::iFun_t iFun=[](Var_t * x,Var_t * v,
            const Var_t * xMin,const Var_t * xMax,const Var_t * vMax) {
        for(size_t i=0;i<N;i++) {
            x->operator[](i)=randD(xMin->operator[](i),xMax->operator[](i));
            v->operator[](i)=0;
        }
    };

    solver_t::fFun_t fFun=[](const Var_t *x,double * fitness) {
        *fitness=10*N;
        auto t=(x->array()*M_2_PI).cos()*10;
        *fitness+=(x->square()-t).sum();
    };

    solver.setPVRange(-5.12,5.12,0.1);

    solver.setiFun(iFun);

    solver.setfFun(fFun);

    solver.setOption(opt);

    solver.initializePop();

    //solver.initialize(iFun,fFun,opt);

    clock_t time=clock();
    solver.run();
    time=clock()-time;

    cout<<"finished in "<<double(time)*1000/CLOCKS_PER_SEC
        <<" miliseconds and "<<solver.generation()<<" generations"<<endl;

    cout<<"result fitness = "<<solver.bestFitness()<<endl;

    const Var_t & result=solver.globalBest().position;
    /*
    cout<<"result = [";
    for(auto i : result) {
        cout<<i<<" , ";
    }
    cout<<"];"<<endl;

    
    cout<<"Trainning Curve=[";
    for(auto i : solver.record()) {
        cout<<i<<" , ";
    }
    cout<<"];"<<endl;

    */

    /*
    cout<<"Population condition:"<<endl;
    for(const auto & i : solver.population()) {
        cout<<"fitness="<<i.fitness<<" , pBest="
           <<i.pBest.fitness<<" , position="<<i.position.transpose()<<" , velocity="<<i.velocity.transpose()<<endl;
    }
    */

}

