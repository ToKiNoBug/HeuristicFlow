#include <iostream>
#include "GABase.h"
#include <ctime>
using namespace std;

void initializeSrand();

int main()
{
    initializeSrand();
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
    return 0;
}

void initializeSrand() {

    std::clock_t c=std::clock();
    if(sizeof(std::clock_t)>sizeof(unsigned int)) {
        std::srand((c>>32)^(c&0xFFFFFFFF));
    } else {
        std::srand(c);
    }
}
