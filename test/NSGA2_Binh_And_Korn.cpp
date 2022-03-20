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


void testNSGA2_Binh_and_Korn() {
    //0<=x_0<=5,  0<=x_1<=3

    using args_t = Heu::BoxNdN<2,Heu::DoubleVectorOption::Std>;

    using solver_t = 
    NSGA2<std::array<double,2>,
            2,
            FITNESS_LESS_BETTER,
            RecordOption::DONT_RECORD_FITNESS,args_t,
            Heu::GADefaults<std::array<double,2>,args_t,Std>::iFunNd,
            nullptr,
            Heu::GADefaults<std::array<double,2>,args_t,Std>::cFunNd<>,
            Heu::GADefaults<std::array<double,2>,args_t,Std>::mFun_d            
            >;

    solver_t algo;

    using Fitness_t = typename solver_t::Fitness_t;

    auto fFun=[](const std::array<double,2> * _x,const args_t *,Fitness_t *f) {
        double f1;
        double f2;
        const double x=_x->operator[](0),y=_x->operator[](1);
        f1=4*(x*x+y*y);
        f2=(x-5)*(x-5)+(y-5)*(y-5);

        double constraint_g1=(x-5)*(x-5)+y*y-25;
        double constraint_g2=7.7-((x-8)*(x-8)+(y+3)*(y+3));

        if(constraint_g1>0) {
            f1=1e4+constraint_g1;
        }

        if(constraint_g2>0) {
            f2=1e4+constraint_g2;
        }

        *f={f1,f2};
    };

    {
        GAOption opt;
        opt.maxGenerations=100;
        opt.populationSize=200;
        opt.maxFailTimes=-1;
        opt.crossoverProb=0.8;
        opt.mutateProb=0.1;
        algo.setOption(opt);
    }
    {
        args_t box;
        box.setMin({0,0});
        box.setMax({5,3});
        box.setLearnRate({0.05,0.03});
        algo.setArgs(box);
    }

    algo.setfFun(fFun);
    algo.initializePop();

    cout<<"Start"<<endl;
    std::clock_t t=std::clock();
    algo.run();
    t=std::clock()-t;
    cout<<"Solving finished in "<<double(t)/CLOCKS_PER_SEC
       <<" seconds and "<<algo.generation()<<"generations"<<endl;

    
    cout<<"paretoFront=[";
    for(const auto & i : algo.pfGenes()) {
        cout<<i->_Fitness[0]<<" , "<<i->_Fitness[1]<<";\n";
    }
    cout<<"];"<<endl;
    
    /*
    cout<<"\n\n\n population=[";
    for(const auto & i : algo.population()) {
        cout<<i.fitness()[0]<<" , "<<i.fitness()[1]<<";\n";
    }
    cout<<"];"<<endl;

    */
    
}
