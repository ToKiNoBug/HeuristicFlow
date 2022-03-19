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


void testTSP(const uint32_t PointNum) {
    static const uint8_t DIM=2;
    //static const double LengthBase=100;
    typedef array<double,DIM> Point_t;
    vector<Point_t> points;
    points.clear();
    points.reserve(PointNum);
    for(uint32_t i=0;i<PointNum;i++) {
        points.emplace_back();
        for(auto & i : points.back()) {
            i=randD();
        }
    }
    /*
    cout<<"Generated random points :"<<endl;
    for(const auto & i : points) {
        cout<<'(';
        for(auto j : i) {
            cout<<j<<" , ";
        }
        cout<<")"<<endl;
    }
    */
    typedef pair<double,uint32_t> permUnit;

    //typedef vector<permUnit> permulation;
    //      var,        less=better,    data src
    SOGA<vector<double>,
            FITNESS_LESS_BETTER,
            DONT_RECORD_FITNESS,
            std::tuple<const vector<Point_t>*>> algo;
    static const uint8_t dataIdx=0;
    typedef tuple<const vector<Point_t>*> Args_t;
    Args_t args;//=make_tuple(PointNum,points.data());
    get<dataIdx>(args)=&points;
    //initialize function

    auto initializeFun=[](vector<double> * x,const Args_t* args) {
        const uint32_t permL=get<dataIdx>(*args)->size();
        x->resize(permL);
        for(uint32_t i=0;i<permL;i++) {
            x->operator[](i)=randD();
        }
    };

    //calculate fitness
    auto calculateFun=[](const vector<double> *x,const Args_t* args,double *f) {

        const uint32_t permL=x->size();
        vector<permUnit> perm(permL);
        for(uint32_t i=0;i<permL;i++) {
            perm[i].first=x->operator[](i);
            perm[i].second=i;
        }
        std::sort(perm.begin(),perm.end(),
                  //compare function for std::sort
                  [](const permUnit &a,const permUnit &b)
                      { return a.first>b.first;});

        double L=0;
        for(uint32_t i=1;i<permL;i++) {
            const Point_t &prev=get<dataIdx>(*args)->operator[](perm[i-1].second);
            const Point_t &cur=get<dataIdx>(*args)->operator[](perm[i].second);
            double curL=0;
            for(uint8_t d=0;d<DIM;d++) {
                curL+=(prev[d]-cur[d])*(prev[d]-cur[d]);
            }
            L+=curL;
        } *f= L;
    };

    //calculate natural perm pathL
    {
        vector<double> natural(PointNum);
        for(uint32_t i=0;i<PointNum;i++) {
            natural[i]=double(i)/PointNum;
        }
        double naturalPathL;
        calculateFun(&natural,&args,&naturalPathL);
        cout<<"default pathL = "<<naturalPathL<<endl;
    }

    auto crossoverFun
            =Heu::GADefaults<vector<double>,Args_t,DoubleVectorOption::Std>::
                cFunXd<Heu::DivCode::Half>;

    auto mutateFun=[](const vector<double>*src,vector<double>*x,const Args_t*) {
        *x=*src;
        const uint32_t permL=x->size();
        if(randD()<0.5) {//modify one element's value
            double & v=x->operator[](std::rand()%permL);
            v+=randD(-0.01,0.01);
            v=std::min(v,1.0);
            v=std::max(v,0.0);
        }
        else {
            const uint32_t flipB=std::rand()%(permL-1);
            const uint32_t flipE=std::rand()%(permL-flipB-1)+flipB+1;
            //cout<<"flipB="<<flipB<<" , flipE"<<flipE<<endl;
            const uint32_t flipN=flipE-flipB+1;
            for(uint32_t flipped=0;flipped*2<=flipN;flipped++) {
                swap(x->operator[](flipB+flipped),x->operator[](flipE-flipped));
            }
        }
    };


    GAOption opt;
    opt.maxGenerations=30*PointNum;
    opt.maxFailTimes=-1;

    algo.setiFun(initializeFun);
    algo.setmFun(mutateFun);
    algo.setfFun(calculateFun);
    algo.setcFun(crossoverFun);
    algo.setOption(opt);
    algo.setArgs(args);
    algo.initializePop();

    cout<<"run!"<<endl;
    std::clock_t c=std::clock();
    algo.run();
    c=std::clock()-c;
    cout<<"finished with "<<algo.generation()<<" generations and "
       <<double(c)/CLOCKS_PER_SEC<<" s\n";
    cout<<"result fitness = "<<algo.bestFitness()<<endl;
}
