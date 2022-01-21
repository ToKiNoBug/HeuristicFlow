#include "def_TestingPSO.h"
#include <cmath>
using namespace OptimT;
using namespace std;

void testPSOBase() {
    PSOBase<std::array<double,10>,
            double,
            RecordOption::DONT_RECORD_FITNESS> * noRec;

    PSOBase<std::array<double,10>,
            double,
            RecordOption::RECORD_FITNESS> * doRec;

    noRec=doRec;

    std::cout<<"sizeof noRec="<<sizeof(typeof (*noRec))<<std::endl;
    std::cout<<"sizeof doRec="<<sizeof(typeof (*doRec))<<std::endl;

    std::array<double,10> a;
    a[0]=0;

    PSO<std::array<double,10>,
            10,
            StdArray,
            FITNESS_GREATER_BETTER,
            DONT_RECORD_FITNESS> psoNoRec;

    PSO<std::array<double,10>,
            10,
            StdArray,
            FITNESS_GREATER_BETTER,
            RECORD_FITNESS> psoDoRec;

    noRec=&psoDoRec;
    doRec=&psoDoRec;

    noRec=&psoNoRec;

    ///doRec=&psoNoRec;
    /// code above failed to compile,
    /// because PSO without fitness recording is not derived from PSOBase with recording
    ///error: incompatible pointer types
    ///assigning to
    /// 'PSOBase<std::array<double, 10>, double, RecordOption::RECORD_FITNESS> *'
    /// from
    /// 'PSO<std::array<double, 10>, 10, StdArray, FITNESS_GREATER_BETTER, DONT_RECORD_FITNESS> *'

}

void testRastriginFun() {
    static const size_t N=10;
    using Var_t = std::array<double,N>;
    using solver_t = PSO<std::array<double,N>,
    N,StdArray,
    FITNESS_LESS_BETTER,
    RECORD_FITNESS>;

    using Args_t = solver_t::Args_t;

    PSOOption opt;
    opt.maxGeneration=1000;
    opt.maxFailTimes=-1;


    solver_t solver;

    solver_t::iFun_t iFun=[](Var_t * x,Var_t * v,
            const Var_t * xMin,const Var_t * xMax,const Var_t * vMax,
            const Args_t *) {
        for(size_t i=0;i<N;i++) {
            x->at(i)=randD(xMin->at(i),xMax->at(i));
            v->at(i)=0;
        }
    };

    solver_t::fFun_t fFun=[](const Var_t *x,const Args_t *,double * fitness) {
        *fitness=10*N;
        for(auto i : *x) {
            //double i=j-1;
            *fitness+=i*i-10*std::cos(M_2_PI*i);
        }
        //cout<<"fFun"<<endl;
    };

    solver.setPVRange(-5.12,5.12,1);

    solver.initialize(iFun,fFun,solver_t::default_ooFun,opt);

    clock_t time=clock();
    solver.run();
    time=clock()-time;

    cout<<"result fitness = "<<solver.bestFitness()<<endl;

    const Var_t & result=solver.globalBest().position;

    cout<<"result = [";
    for(auto i : result) {
        cout<<i<<" , ";
    }
    cout<<"];"<<endl;


}
