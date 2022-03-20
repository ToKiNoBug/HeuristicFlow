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

///Zitzler–Deb–Thiele's function N. 3
void testNSGA2_ZDT3() {
    //0<=x_i<=1, 1<=i<=30
    static const size_t XNum=30;
    static const double r=0.2;

    Heu::NSGA2<std::array<double,XNum>,
            2,
            FITNESS_LESS_BETTER,
            RECORD_FITNESS> algo;

    void (*iFun)(std::array<double,XNum>*) =
    [] (std::array<double,XNum> * x) {
        for(size_t i=0;i<XNum;i++) {
            x->operator[](i)=Heu::randD();
        }
    };

    void (*fFun)
            (const std::array<double,XNum>* x, Eigen::Array<double,2,1> *f)
            =[](const std::array<double,XNum>* x, Eigen::Array<double,2,1> *f) {
      f->operator[](0)=x->operator[](0);
      const double && f1=std::move(f->operator[](0));
      double g=0;
      for(size_t i=1;i<XNum;i++) {
          g+=x->operator[](i);
      }
      g=1+g*9.0/(XNum-1);
      double && f1_div_g=f1/g;
      double && h=1-std::sqrt(f1_div_g)-f1_div_g*std::sin(10*M_PI*f1);
      f->operator[](1)=g*h;
    };

    auto cFun=
            Heu::GADefaults<std::array<double,XNum>>::
                cFunNd<(Heu::DivEncode<1,2>::code)>;

    void (*mFun)(const std::array<double,XNum>*src,std::array<double,XNum>*)=
            [](const std::array<double,XNum>*src,std::array<double,XNum>*x){
        *x=*src;
        const size_t mutateIdx=Heu::randIdx(XNum);

        x->operator[](mutateIdx)+=0.005*Heu::randD(-1,1);

        x->operator[](mutateIdx)=std::min(x->operator[](mutateIdx),1.0);
        x->operator[](mutateIdx)=std::max(x->operator[](mutateIdx),0.0);
    };

    GAOption opt;
    opt.populationSize=2000;
    opt.maxFailTimes=500;
    opt.maxGenerations=2000;

    algo.setiFun(iFun);
    algo.setmFun(mFun);
    algo.setfFun(fFun);
    algo.setcFun(cFun);
    algo.setOption(opt);
    algo.initializePop();

    cout<<"Start"<<endl;
    std::clock_t t=std::clock();
    algo.run();
    t=std::clock()-t;
    cout<<"Solving finished in "<<double(t)/CLOCKS_PER_SEC
       <<"seconds and "<<algo.generation()<<"generations"<<endl;
    std::vector<Eigen::Array<double,2,1>> paretoFront;
    algo.paretoFront(paretoFront);
    cout<<"paretoFront=[";
    for(const auto & i : paretoFront) {
        cout<<i[0]<<" , "<<i[1]<<";\n";
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

void testNSGA2_Kursawe() {
    //1<=i<=3,  -5<=x_i<=5
    NSGA2<std::array<double,3>,
            2,
            FITNESS_LESS_BETTER,
            DONT_RECORD_FITNESS> algo;
    auto iFun=[](std::array<double,3> * x) {
        for(auto & i : *x) {
            i=randD(-5,5);
        }
    };
    auto fFun=[](const std::array<double,3> * x,Eigen::Array<double,2,1> *f) {
        double f1=0,f2=0;
        for(int i=0;i<2;i++) {
            f1+=-10*exp(-0.2*sqrt((x->operator[](i))*(x->operator[](i))
                    +(x->operator[](i+1))*(x->operator[](i+1))));
        }
        for(int i=0;i<3;i++) {
            f2+=pow(abs(x->operator[](i)),0.8)+5*sin(x->operator[](i)*x->operator[](i)*x->operator[](i));
        }
        f->operator[](0)=f1;
        f->operator[](1)=f2;
    };

    auto cFun=[](const std::array<double,3> *p1,const std::array<double,3> *p2,
            std::array<double,3> *ch1,std::array<double,3> *ch2) {
        for(int i=0;i<3;i++) {
            static const double r=0.2;
            ch1->operator[](i)=r*p1->operator[](i)+(1-r)*p2->operator[](i);
            ch2->operator[](i)=r*p2->operator[](i)+(1-r)*p1->operator[](i);
        }
    };

    auto mFun=[](const std::array<double,3> * src,std::array<double,3> * x) {
        *x=*src;
        const size_t idx=randIdx(3);
        x->operator[](idx)+=0.1*randD(-1,1);
        x->operator[](idx)=std::min(x->operator[](idx),5.0);
        x->operator[](idx)=std::max(x->operator[](idx),-5.0);
    };

    GAOption opt;
    opt.maxGenerations=200;
    opt.populationSize=100;
    opt.maxFailTimes=-1;

    algo.setiFun(iFun);
    algo.setmFun(mFun);
    algo.setfFun(fFun);
    algo.setcFun(cFun);
    algo.setOption(opt);
    algo.initializePop();

    cout<<"Start"<<endl;
    std::clock_t t=std::clock();
    algo.run();
    t=std::clock()-t;
    cout<<"Solving finished in "<<double(t)/CLOCKS_PER_SEC
       <<" seconds and "<<algo.generation()<<" generations"<<endl;
    std::vector<Eigen::Array<double,2,1>> paretoFront;
    algo.paretoFront(paretoFront);
    cout<<"paretoFront=[";
    for(const auto & i : paretoFront) {
        cout<<i[0]<<" , "<<i[1]<<";\n";
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
