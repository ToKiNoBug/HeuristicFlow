/*
 Copyright © 2021  TokiNoBug
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

#include "def_TestingFuns.h"
#include <Genetic>
#include <iostream>
#include <Eigen/Dense>
using namespace OptimT;
using namespace std;

void testAckley_withRecord() {
    GAOption opt;
    opt.maxGenerations=1000;
    opt.maxFailTimes=100;
    static const uint8_t MinIdx=0,MaxIdx=1,LrIdx=2;
    SOGA<array<double,2>,false,true,array<double,2>,array<double,2>,double> algo;

    /*
    GABase<array<double,2>,double,true,array<double,2>,array<double,2>,double> * basePtr;
    basePtr=&algo;
    */

    tuple<array<double,2>,array<double,2>,double> args;
    get<MinIdx>(args)={-5,-5};
    get<MaxIdx>(args)={5,5};
    get<LrIdx>(args)=0.05;

    algo.initialize(
                //initializer for array<double,2>
    [](array<double,2>* x,const typeof(args) * a) {
        for(uint32_t idx=0;idx<x->size();idx++) {
            x->operator[](idx)=randD(get<MinIdx>(*a)[idx],get<MaxIdx>(*a)[idx]);
        }},
    //Ackely function
    [](const array<double,2>* _x,const typeof(args) *) {
        double x=_x->operator[](0),y=_x->operator[](1);
        return -20*exp(-0.2*sqrt(0.5*(x*x+y*y)))
                -exp(0.5*(cos(M_2_PI*x)+cos(M_2_PI*y)))
                +20+M_E;},
    //crossover
    [](const array<double,2>* x,const array<double,2>* y,
            array<double,2> *X,array<double,2>*Y,
            const typeof(args) *) {
        const array<double,2> &copyx=*x,&copyy=*y;
        for(uint32_t idx=0;idx<x->size();idx++) {
            if(rand()%2)
                X->operator[](idx)=copyy[idx];
            else {
                X->operator[](idx)=copyx[idx];
            }
            if(rand()%2)
                Y->operator[](idx)=copyx[idx];
            else {
                Y->operator[](idx)=copyy[idx];
            }
        }},
    //mutate
    [](array<double,2>* x,const typeof(args) * a) {
        uint8_t mutateIdx=rand()%x->size();
        x->operator[](mutateIdx)+=get<LrIdx>(*a)*randD(-1,1);
        if(rand()%2)
            x->operator[](mutateIdx)*=-1;
        for(uint32_t idx=0;idx<x->size();idx++) {
            auto & var=x->operator[](idx);
            var=min(var,get<MaxIdx>(*a)[idx]);
            var=max(var,get<MinIdx>(*a)[idx]);
        }
    },
    //no other options
    nullptr,
    opt,
    args);

    algo.run();

    cout<<"Solving spend "<<algo.generation()<<" generations\n";
    cout<<"Result = ["<<algo.result()[0]<<" , "<<algo.result()[1]<<"]\n";

    cout<<"Fitness history :\n";

    for(auto i : algo.record()) {
        cout<<i<<'\n';
    }
    cout<<endl;

    //cout<<"Result idx = "<<algo.eliteIdx()<<endl;
}


void testSingleNumber() {
    SOGA<double,true,false> algo;
    GAOption opt;
    opt.crossoverProb=0.8;
    opt.mutateProb=0.05;
    opt.maxFailTimes=100;
    opt.populationSize=100;
    opt.maxGenerations=3000;
    algo.initialize(
        [](double * x,const tuple<>*){*x=randD();},
        [](const double * x,const tuple<>*){return 1.0/(abs(*x-3.0)+1e-8);},
        [](const double * x,const double *y,double *X,double *Y,const tuple<>*)
            {double mean=(*x+*y)/2;*X=(*x+mean)/2;*Y=(*y+mean)/2;},
        [](double * x,const tuple<>*)
            {
        //if(rand()%2)*x*=-1;
             *x+=randD(-1,1)*0.5;},
    nullptr,
    opt,
    tuple<>()
    );

    cout << "Start" << endl;
    algo.run();
    cout<<"Finished with "<<algo.generation()<<" generations"<<endl;
    cout<<"result = "<<algo.result()<<endl;
}

void testWithEigenLib() {
    Eigen::Array4d target(1,2,3,4);
    target/=target.sum();
    cout<<"Target="<<target.transpose()<<endl;
    GAOption opt;
    opt.populationSize=100;
    opt.maxGenerations=3000;
    opt.maxFailTimes=50;
    opt.crossoverProb=0.8;
    opt.mutateProb=0.05;

    SOGA<Eigen::Array4d,true,false,Eigen::Array4d,Eigen::Array4d,Eigen::Array4d,double> algo;
    //val min max learning_rate
    tuple<Eigen::Array4d,Eigen::Array4d,Eigen::Array4d,double> Arg;
    static const uint8_t TargetOffset=0,MinOffset=1,MaxOffset=2,LROffset=3;
    get<TargetOffset>(Arg)=target;
    get<MinOffset>(Arg).setConstant(-1);
    get<MaxOffset>(Arg).setConstant(2);
    get<LROffset>(Arg)=0.01;

    //typeof(Arg)


    algo.initialize(
                    [](Eigen::Array4d* x,const typeof(Arg)*){x->setRandom();},
    [](const Eigen::Array4d* x,const typeof(Arg)* arg){
        return -(*x-get<TargetOffset>(*arg)).square().maxCoeff();
    },
    [](const Eigen::Array4d*x,const Eigen::Array4d*y,
            Eigen::Array4d*X,Eigen::Array4d*Y,
            const typeof(Arg)*) {
        for(uint32_t i=0;i<4;i++) {
            X->operator()(i)=
                    (rand()%2)?
                        x->operator()(i):y->operator()(i);
            Y->operator()(i)=
                    (rand()%2)?
                        x->operator()(i):y->operator()(i);
        }},
    [](Eigen::Array4d*x,const typeof(Arg)* arg) {
        uint32_t idx=rand()%4;
        x->operator()(idx)+=randD(-1,1)*get<LROffset>(*arg);

        if(x->operator()(idx)>get<MaxOffset>(*arg)(idx)) {
            x->operator()(idx)=get<MaxOffset>(*arg)(idx);
        }
        if(x->operator()(idx)<get<MinOffset>(*arg)(idx)) {
            x->operator()(idx)=get<MinOffset>(*arg)(idx);
        }},
    nullptr,
    opt,
    Arg);


    algo.run();
    cout<<"Solving spend "<<algo.generation()<<" generations\n";
    cout<<"Result = "<<algo.result().transpose()<<endl;
}

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
    cout<<"Generated random points :"<<endl;
    for(const auto & i : points) {
        cout<<'(';
        for(auto j : i) {
            cout<<j<<" , ";
        }
        cout<<")"<<endl;
    }
    typedef pair<double,uint32_t> permUnit;

    //typedef vector<permUnit> permulation;
    //      var,        less=better,    data src
    SOGA<vector<double>,false,false,const vector<Point_t>*> algo;
    static const uint8_t dataIdx=0;
    typedef tuple<const vector<Point_t>*> Args_t;
    Args_t args;//=make_tuple(PointNum,points.data());
    get<dataIdx>(args)=&points;
    //initialize function

    auto initializeFun=[](vector<double> * x,const Args_t* args) {
        const uint32_t permL=get<dataIdx>(*args)->size();
        x->resize(permL);
        for(uint32_t i=0;i<permL;i++) {
            x->at(i)=OptimT::randD();
        }
    };

    //calculate fitness
    auto calculateFun=[](const vector<double> *x,const Args_t* args) {
        const uint32_t permL=x->size();
        vector<permUnit> perm(permL);
        for(uint32_t i=0;i<permL;i++) {
            perm[i].first=x->at(i);
            perm[i].second=i;
        }
        std::sort(perm.begin(),perm.end(),
                  //compare function for std::sort
                  [](const permUnit &a,const permUnit &b)
                      { return a.first>b.first;});
        double L=0;
        for(uint32_t i=1;i<permL;i++) {
            const Point_t &prev=get<dataIdx>(*args)->at(perm[i-1].second);
            const Point_t &cur=get<dataIdx>(*args)->at(perm[i].second);
            double curL=0;
            for(uint8_t d=0;d<DIM;d++) {
                curL+=(prev[d]-cur[d])*(prev[d]-cur[d]);
            }
            L+=curL;
        } return L;};

    //calculate natural perm pathL
    {
        vector<double> natural(PointNum);
        for(uint32_t i=0;i<PointNum;i++) {
            natural[i]=double(i)/PointNum;
        }
        double naturalPathL=calculateFun(&natural,&args);
        cout<<"default pathL = "<<naturalPathL<<endl;
    }

    auto crossoverFun=[](const vector<double>*p1,const vector<double>*p2,
                      vector<double>*c1,vector<double>*c2,const Args_t*) {
        const uint32_t permL=p1->size();
        c1->resize(permL);c2->resize(permL);
        for(uint32_t i=0;i<permL;i++) {
            c1->at(i)=(std::rand()%2)?p1->at(i):p2->at(i);
            c2->at(i)=(std::rand()%2)?p1->at(i):p2->at(i);
        }
    };

    auto mutateFun=[](vector<double>*x,const Args_t*) {
        const uint32_t permL=x->size();
        if(randD()<0.5) {//modify one element's value
            double & v=x->at(std::rand()%permL);
            v+=randD(-0.01,0.01);
            v=min(v,1.0);
            v=max(v,0.0);
        }
        else {
            const uint32_t flipB=std::rand()%(permL-1);
            const uint32_t flipE=std::rand()%(permL-flipB-1)+flipB+1;
            //cout<<"flipB="<<flipB<<" , flipE"<<flipE<<endl;
            const uint32_t flipN=flipE-flipB+1;
            for(uint32_t flipped=0;flipped*2<=flipN;flipped++) {
                swap(x->at(flipB+flipped),x->at(flipE-flipped));
            }
        }
    };


    GAOption opt;
    opt.maxGenerations=30*PointNum;
    algo.initialize(initializeFun,calculateFun,crossoverFun,mutateFun,nullptr,
                    opt,args);

    cout<<"run!"<<endl;
    algo.run();
    cout<<"finished with "<<algo.generation()<<" generations\n";
    cout<<"result fitness = "<<algo.eliteIt()->fitness()<<endl;



}

