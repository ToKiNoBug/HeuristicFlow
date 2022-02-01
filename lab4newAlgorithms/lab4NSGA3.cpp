#include "lab4NSGA3.h"
#include <iostream>
using namespace std;


void pri_makeRP(const int64_t dimN,const int64_t precision,
    const int64_t curDim,const int64_t curP,const int64_t accum,
    Eigen::ArrayXd & rec,
    vector<Eigen::ArrayXd> & dst) {
    
    if(curDim+1>=dimN) {
        rec[dimN-1]=1.0-double(accum)/precision;
        dst.push_back(rec);
        return;
    }
    

    for(int64_t p=0;p+accum<=precision;p++) {
        if(curDim>=0)
            rec[curDim]=double(p)/precision;
        pri_makeRP(dimN,precision,curDim+1,p,accum+p,rec,dst);
    }
}

void pri_startRP(const int64_t dimN,const int64_t precision,vector<Eigen::ArrayXd> & dst) {
    Eigen::ArrayXd rec;
    rec.setConstant(dimN,1,4);
    pri_makeRP(dimN,precision,0,0,0,rec,dst);
}

vector<Eigen::ArrayXd> makeReferencePoints(const int64_t dimN,const int64_t precision) {
    if(precision<=0) {
        exit(114514);
    }
    vector<Eigen::ArrayXd> points;

    pri_startRP(dimN,precision,points);

    cout<<"RPS=[\n";
    for(const auto & i : points) {
        cout<<double(precision)*i.transpose()<<'\n';
    }
    cout<<"];"<<endl;

    return points;

}
