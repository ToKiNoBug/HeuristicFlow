/*
 Copyright © 2021  TokiNoBug
This file is part of AlgoTemplates.

    AlgoTemplates is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AlgoTemplates is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AlgoTemplates.  If not, see <https://www.gnu.org/licenses/>.

*/

#include <iostream>
#include "Genetic.h"
#include <ctime>

#include <Eigen/Dense>

using namespace std;

using namespace AT;

void initializeSrand();

void testSingleNumber();

void testWithEigen();

void testAckley();

int main()
{
    initializeSrand();
    //testSingleNumber();
    testAckley();
    return 0;
}

void initializeSrand() {
    std::time_t t=std::time(nullptr);
    std::srand((t>>32)^(t&0xFFFFFFFF));
    /*
    std::clock_t c=std::clock();
    if(sizeof(std::clock_t)>sizeof(unsigned int)) {
        std::srand((c>>32)^(c&0xFFFFFFFF));
    } else {
        std::srand(c);
    }*/
}

void testSingleNumber() {
    GA<double,true> algo;
    GAOption opt;
    opt.crossoverProb=0.8;
    opt.mutateProb=0.05;
    opt.maxFailTimes=100;
    opt.populationSize=100;
    opt.maxGenerations=3000;
    algo.initialize(
        [](double * x,const tuple<>*){*x=randD();},
        [](const double * x,const tuple<>*){return 1.0/(std::abs(*x-3.0)+1e-8);},
        [](const double * x,const double *y,double *X,double *Y,const tuple<>*)
            {double mean=(*x+*y)/2;*X=(*x+mean)/2;*Y=(*y+mean)/2;},
        [](double * x,const tuple<>*)
            {
        //if(std::rand()%2)*x*=-1;
             *x+=randD(-1,1)*0.5;},
    nullptr,
    opt,
    std::tuple<>()
    );

    cout << "Start" << endl;
    algo.run();
    cout<<"Finished with "<<algo.generation()<<" generations"<<endl;
    cout<<"result = "<<algo.result()<<endl;
}
/*
void testWithEigen() {
    Eigen::Array4d target(1,2,3,4);
    target/=target.sum();
    cout<<"Target="<<target.transpose()<<endl;
    GAOption opt;
    opt.populationSize=100;
    opt.maxGenerations=3000;
    opt.maxFailTimes=50;
    opt.crossoverProb=0.8;
    opt.mutateProb=0.05;

    GA<Eigen::Array4d,true,Eigen::Array4d,Eigen::Array4d,Eigen::Array4d,double> algo;
    //val min max learning_rate
    std::tuple<Eigen::Array4d,Eigen::Array4d,Eigen::Array4d,double> Arg;
    static const uint8_t TargetOffset=0,MinOffset=1,MaxOffset=2,LROffset=3;
    std::get<TargetOffset>(Arg)=target;
    std::get<MinOffset>(Arg).setConstant(-1);
    std::get<MaxOffset>(Arg).setConstant(2);
    std::get<LROffset>(Arg)=0.01;

    //typeof(Arg)


    algo.initialize(
                    [](Eigen::Array4d* x,const typeof(Arg)*){x->setRandom();},
    [](const Eigen::Array4d* x,const typeof(Arg)* arg){
        return -(*x-std::get<TargetOffset>(*arg)).square().maxCoeff();
    },
    [](Eigen::Array4d*x,Eigen::Array4d*y,const typeof(Arg)*) {
        for(uint32_t i=0;i<4;i++) {
            if(std::rand()%2==0) {
                std::swap(x->operator()(i),y->operator()(i));
            }
        }},
    [](Eigen::Array4d*x,const typeof(Arg)* arg) {
        uint32_t idx=std::rand()%4;
        x->operator()(idx)+=randD(-1,1)*std::get<LROffset>(*arg);

        if(x->operator()(idx)>std::get<MaxOffset>(*arg)(idx)) {
            x->operator()(idx)=std::get<MaxOffset>(*arg)(idx);
        }
        if(x->operator()(idx)<std::get<MinOffset>(*arg)(idx)) {
            x->operator()(idx)=std::get<MinOffset>(*arg)(idx);
        }},
    nullptr,
    opt,
    Arg);


    algo.run();
    cout<<"Solving spend "<<algo.generation()<<" generations\n";
    cout<<"Result = "<<algo.result().transpose()<<endl;
    cout<<"Result idx = "<<algo.eliteIdx()<<endl;
}
*/

void testAckley() {
    GAOption opt;
    opt.maxGenerations=1000;
    opt.maxFailTimes=100;
    static const uint8_t MinIdx=0,MaxIdx=1,LrIdx=2;
    GA<std::array<double,2>,false,std::array<double,2>,std::array<double,2>,double> algo;

    std::tuple<std::array<double,2>,std::array<double,2>,double> args;
    std::get<MinIdx>(args)={-5,-5};
    std::get<MaxIdx>(args)={5,5};
    std::get<LrIdx>(args)=0.05;

    algo.initialize(
                //initializer for std::array<double,2>
    [](std::array<double,2>* x,const typeof(args) * a) {
        for(uint32_t idx=0;idx<x->size();idx++) {
            x->operator[](idx)=randD(std::get<MinIdx>(*a)[idx],std::get<MaxIdx>(*a)[idx]);
        }},
    //Ackely function
    [](const std::array<double,2>* _x,const typeof(args) *) {
        double x=_x->operator[](0),y=_x->operator[](1);
        return -20*exp(-0.2*sqrt(0.5*(x*x+y*y)))
                -exp(0.5*(cos(M_2_PI*x)+cos(M_2_PI*y)))
                +20+M_E;},
    //crossover
    [](const std::array<double,2>* x,const std::array<double,2>* y,
            std::array<double,2> *X,std::array<double,2>*Y,
            const typeof(args) *) {
        const std::array<double,2> &copyx=*x,&copyy=*y;
        for(uint32_t idx=0;idx<x->size();idx++) {
            if(std::rand()%2)
                X->operator[](idx)=copyy[idx];
            else {
                X->operator[](idx)=copyx[idx];
            }
            if(std::rand()%2)
                Y->operator[](idx)=copyx[idx];
            else {
                Y->operator[](idx)=copyy[idx];
            }
        }},
    //mutate
    [](std::array<double,2>* x,const typeof(args) * a) {
        uint8_t mutateIdx=std::rand()%x->size();
        x->operator[](mutateIdx)+=std::get<LrIdx>(*a)*randD(-1,1);
        if(std::rand()%2)
            x->operator[](mutateIdx)*=-1;
        for(uint32_t idx=0;idx<x->size();idx++) {
            auto & var=x->operator[](idx);
            var=std::min(var,std::get<MaxIdx>(*a)[idx]);
            var=std::max(var,std::get<MinIdx>(*a)[idx]);
        }
    },
    //no other options
    nullptr,
    opt,
    args);

    algo.run();

    cout<<"Solving spend "<<algo.generation()<<" generations\n";
    cout<<"Result = ["<<algo.result()[0]<<" , "<<algo.result()[1]<<"]\n";
    //cout<<"Result idx = "<<algo.eliteIdx()<<endl;
}

