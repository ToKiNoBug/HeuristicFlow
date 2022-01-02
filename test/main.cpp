#include <iostream>
#include "GABase.h"
#include <ctime>

#include <Eigen/Dense>

using namespace std;

void initializeSrand();

void testSingleNumber();

void testWithEigen();

int main()
{
    initializeSrand();
    testWithEigen();
    return 0;
}

void initializeSrand() {

    std::srand(std::clock());
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
    opt.maxFailTimes=50;
    opt.populationSize=100;
    opt.maxGenerations=300;
    algo.initialize(opt,
        [](double * x,const tuple<>*){*x=randD();},
        [](const double * x,const tuple<>*){return 1.0/(std::abs(*x-3.0)+1e-8);},
        [](double * x,double *y,const tuple<>*)
            {double mean=(*x+*y)/2;*x=(*x+mean)/2;*y=(*y+mean)/2;},
        [](double * x,const tuple<>*)
            {if(std::rand()%2)*x*=-1;
             *x+=randD()*0.1;},
    std::tuple<>()
    );

    cout << "Start" << endl;
    algo.run();
    cout<<"Finished"<<endl;
    cout<<"result = "<<algo.result()<<endl;
    cout<<"result idx = "<<algo.eliteIdx()<<endl;
}

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


    algo.initialize(opt,
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
        /*
        if(std::rand()%2==0) {
            x->operator()(idx)*=-1;
        }*/
        if(x->operator()(idx)>std::get<MaxOffset>(*arg)(idx)) {
            x->operator()(idx)=std::get<MaxOffset>(*arg)(idx);
        }
        if(x->operator()(idx)<std::get<MinOffset>(*arg)(idx)) {
            x->operator()(idx)=std::get<MinOffset>(*arg)(idx);
        }
    },
    Arg);


    algo.run();

    cout<<"Result = "<<algo.result().transpose()<<endl;
    cout<<"Result idx = "<<algo.eliteIdx()<<endl;
}
