#include "def_TestingPSO.h"
#include <cmath>
#include <iostream>
using namespace OptimT;
using namespace std;

void testPSOBase() {
    using Var_t = std::vector<double>;
    PSOAbstract<Var_t,
            double,
            RecordOption::DONT_RECORD_FITNESS> * abstractNoRec;

    PSOAbstract<Var_t,
            double,
            RecordOption::RECORD_FITNESS> * abstractDoRec;

    PSOBase<Var_t,10,
            double,
            RecordOption::DONT_RECORD_FITNESS> * baseNoRec;

    PSOBase<Var_t,0,
            double,
            RecordOption::RECORD_FITNESS> * baseDoRec;

    //DoRec is derived from NoRec
    abstractNoRec=abstractDoRec;

    //base is derived from abstract
    abstractNoRec=baseNoRec;
    abstractDoRec=baseDoRec;




    



}

void testRastriginFun() {
    static const size_t N=50;
    using solver_t = PSO_Eigen<N,
    FITNESS_LESS_BETTER,
    RECORD_FITNESS>;
    
    using Var_t = Eigen::Array<double,N,1>;

    using Args_t = solver_t::Args_t;

    PSOOption opt;
    opt.maxGeneration=10000;
    opt.maxFailTimes=opt.maxGeneration/10;
    opt.inertiaFactor=0.8;
    opt.learnFactorG=2;
    opt.learnFactorP=2;


    solver_t solver;

    solver_t::iFun_t iFun=[](Var_t * x,Var_t * v,
            const Var_t * xMin,const Var_t * xMax,const Var_t * vMax,
            const Args_t *) {
        for(size_t i=0;i<N;i++) {
            x->operator[](i)=randD(xMin->operator[](i),xMax->operator[](i));
            v->operator[](i)=0;
        }
    };

    solver_t::fFun_t fFun=[](const Var_t *x,const Args_t *,double * fitness) {
        *fitness=10*N;
        *fitness+=(x->square()-10*(*x*M_2_PI).cos()).sum();

    };

    solver.setPVRange(-5.12,5.12,1);

    solver.initialize(iFun,fFun,solver_t::default_ooFun,opt);

    clock_t time=clock();
    solver.run();
    time=clock()-time;

    cout<<"finished in "<<double(time)*1000/CLOCKS_PER_SEC
        <<" miliseconds and "<<solver.generation()<<" generations"<<endl;

    cout<<"result fitness = "<<solver.bestFitness()<<endl;

    const Var_t & result=solver.globalBest().position;

    cout<<"result = [";
    for(auto i : result) {
        cout<<i<<" , ";
    }
    cout<<"];"<<endl;


}

