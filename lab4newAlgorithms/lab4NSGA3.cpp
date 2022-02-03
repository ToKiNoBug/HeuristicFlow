#include "lab4NSGA3.h"
#include <iostream>
#include <OptimTemplates/Global>

using namespace std;


void pri_makeRP(const uint64_t dimN,const uint64_t precision,
    const uint64_t curDim,const uint64_t curP,const uint64_t accum,
    Eigen::ArrayXd & rec,
    vector<Eigen::ArrayXd> & dst) {
    
    if(curDim+1>=dimN) {
        rec[dimN-1]=1.0-double(accum)/precision;
        dst.push_back(rec);
        return;
    }
    

    for(uint64_t p=0;p+accum<=precision;p++) {
        if(curDim>=0)
            rec[curDim]=double(p)/precision;
        pri_makeRP(dimN,precision,curDim+1,p,accum+p,rec,dst);
    }
}

void pri_startRP(const uint64_t dimN,const uint64_t precision,vector<Eigen::ArrayXd> & dst) {
    Eigen::ArrayXd rec;
    rec.setConstant(dimN,1,4);
    pri_makeRP(dimN,precision,0,0,0,rec,dst);
}

vector<Eigen::ArrayXd> makeReferencePoints(const uint64_t dimN,const uint64_t precision) {
    if(precision<=0) {
        exit(114514);
    }
    vector<Eigen::ArrayXd> points;

    points.reserve(OptimT::NchooseK(dimN+precision-1,precision));

    pri_startRP(dimN,precision,points);

    /*
    cout<<"RPS=[\n";
    for(const auto & i : points) {
        cout<<i.transpose()<<'\n';
    }
    cout<<"];"<<endl;
    */

    return points;

}

Eigen::ArrayXd sample2Intercept(Eigen::MatrixXd P) {
    auto P_transpose_inv=P.transpose().inverse();
    auto ONE=Eigen::Matrix<double,Eigen::Dynamic,1>::Ones(P.cols(),1);
    auto one_div_intercept=(P_transpose_inv*ONE).array();
    return 1.0/one_div_intercept;
}

void testNSGA3Expri() {
    testNSGA3 solver;

    solver.setPrecision(4);
    OptimT::GAOption opt;
    opt.maxGenerations=1000;
    opt.maxFailTimes=400;
    opt.populationSize=400;
    opt.crossoverProb=0.8;
    opt.mutateProb=0.05;

    solver.initialize(testNSGA3::iFun,testNSGA3::fFun,testNSGA3::cFun,testNSGA3::mFun,nullptr,opt);
    clock_t c=clock();
    solver.run();
    c=clock()-c;

    cout<<"Solving finished in "<<c<<" ms with "<<solver.generation()<<" generations"<<endl;

    vector<Eigen::Array3d> PF;
    solver.paretoFront(PF);

    cout<<"PFV=[\n";
    for(auto & i : PF) {
        cout<<i.transpose()<<'\n';
    }
    cout<<"];\n"<<endl;

}