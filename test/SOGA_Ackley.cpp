// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Eigen/Dense>
#include <HeuristicFlow/Genetic>
#include <iostream>
#include <ctime>
using namespace Heu;
using namespace std;

void testAckley_withRecord() {

    using args_t = Heu::BoxNdS<2,Std>;

    using solver_t = 
    SOGA<array<double,2>,
            Heu::FITNESS_LESS_BETTER,
            Heu::RECORD_FITNESS,
            args_t,
    Heu::GADefaults<array<double,2>,args_t,Std>::iFunNd<>,
    nullptr,
    Heu::GADefaults<array<double,2>,args_t,Std>::cFunNd,
    Heu::GADefaults<array<double,2>,args_t,Std>::mFun_d<>>;
    solver_t algo;
    
    GAOption opt;
    opt.populationSize=50;
    opt.maxFailTimes=-1;
    opt.maxGenerations=100;
    
    algo.setOption(opt);

    {
        args_t args;
        args.setMin(-5);
        args.setMax(5);
        args.setLearnRate(0.05);
        algo.setArgs(args);
    }

    algo.setfFun(
    //Ackely function
    [](const array<double,2>* _x,const solver_t::ArgsType *,double * f) {
        double x=_x->operator[](0),y=_x->operator[](1);
        *f= -20*exp(-0.2*sqrt(0.5*(x*x+y*y)))
                -exp(0.5*(cos(M_2_PI*x)+cos(M_2_PI*y)))
                +20+M_E;}
    );

    algo.initializePop();

    std::clock_t t=std::clock();
    algo.run();
    t=std::clock()-t;
    //cout<<algo.bestFitness();

    cout<<"Solving spend "<<algo.generation()<<" generations in "
       <<double(t)/CLOCKS_PER_SEC<<" sec\n";
    cout<<"Result = ["<<algo.result()[0]<<" , "<<algo.result()[1]<<"]\n";

    
    cout<<"Fitness history :\n";

    for(auto i : algo.record()) {
        cout<<i<<'\n';
    }
    cout<<endl;
    
}

int main()
{
    testAckley_withRecord();
    system("pause");
    return 0;
}