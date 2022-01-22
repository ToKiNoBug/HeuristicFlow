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
    static const size_t N=10;
    using solver_t = PSO_Eigen<N,
    FITNESS_LESS_BETTER,
    RECORD_FITNESS>;
    
    using Var_t = EigenVar_t<N>;;

    using Args_t = solver_t::Args_t;

    PSOOption opt;
    opt.maxGeneration=50*N;
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

    solver.setPVRange(-5.12,5.12,0.1);

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

    
    cout<<"Trainning Curve=[";
    for(auto i : solver.record()) {
        cout<<i<<" , ";
    }
    cout<<"];"<<endl;

    /*
    cout<<"Population condition:"<<endl;
    for(const auto & i : solver.population()) {
        cout<<"fitness="<<i.fitness<<" , pBest="
           <<i.pBest.fitness<<" , position="<<i.position.transpose()<<" , velocity="<<i.velocity.transpose()<<endl;
    }
    */

}

#include <algorithm>

void testTSP(const size_t N) {
static const size_t SpaceDim=2;

    using Var_t = Eigen::ArrayXd;

    using DistanceMat_t = Eigen::ArrayXXd;

    using Solver_t = PSO_Eigen<Dynamic,
    FITNESS_LESS_BETTER,
    RECORD_FITNESS,DistanceMat_t>;


    cout<<Enum2String(Solver_t::Flag)<<endl;

    using Args_t = Solver_t::Args_t;

    Solver_t solver;

    using sortUnit = std::pair<double,size_t>;

    Solver_t::fFun_t fFun=[](const Var_t* x,const Args_t * args,double * fitness) {
        const size_t N=std::get<0>(*args).rows();
        std::vector<sortUnit> sortSpace(N);
        for(size_t i=0;i<N;i++) {
            sortSpace[i]=std::make_pair(x->operator[](i),i);
        }

        
        static const auto cmpFun=[](const sortUnit & a,const sortUnit & b) {
            return a.first<b.first;
        };


        std::sort(sortSpace.begin(),sortSpace.end(),cmpFun);

        *fitness=0;

        for(size_t i=0;i+1<N;i++) {
            *fitness+=std::get<0>(*args)(sortSpace[i].second,sortSpace[i+1].second);
        }
    };

    DistanceMat_t dMat(N,N);
    
    dMat.setZero();

    {
        Eigen::Array<double,SpaceDim,Eigen::Dynamic> points;
        points.setRandom(SpaceDim,N);

        for(size_t r=0;r<N;r++) {
            for(size_t c=0;c<N;c++) {
                if(r==c) continue;
                dMat(r,c)=(points.col(r)-points.col(c)).square().sum();
            }
        }
    }

    PSOOption opt;
    opt.inertiaFactor=0.8;
    opt.learnFactorG=2;
    opt.learnFactorP=2;
    opt.maxGeneration=50*N;
    opt.maxFailTimes=opt.maxGeneration/10;
    opt.populationSize=100;

    

    solver.setDimensions(N);

    solver.setPVRange(0,1,0.5);

    solver.initialize(Solver_t::default_iFun,
                        fFun,
                        Solver_t::default_ooFun,
                        opt,
                        std::make_tuple(dMat));

    {
        double iniFitness;
        fFun(&solver.population().front().position,&solver.args(),&iniFitness);
        cout<<"iniFitness = "<<iniFitness<<endl;
    }

    clock_t c=clock();
    solver.run();
    c=clock()-c;

    cout<<"finished in "<<double(c)*1000/CLOCKS_PER_SEC
        <<" miliseconds and "<<solver.generation()<<" generations"<<endl;

    cout<<"result fitness = "<<solver.bestFitness()<<endl;

    cout<<"Trainning Curve=[";
    for(auto i : solver.record()) {
        cout<<i<<" , ";
    }
    cout<<"];"<<endl;
    
    cout<<"Population condition:"<<endl;
    for(const auto & i : solver.population()) {
        cout<<"fitness="<<i.fitness<<" , pBest="<<i.pBest.fitness<<endl;
    }
    



}
