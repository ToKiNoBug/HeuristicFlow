#include "def_TestingMatrix.h"
#include <OptimTemplates/SimpleMatrix>
#include <stdint.h>
#include <iostream>
#include <Eigen/Dense>
#include <array>
using namespace OptimT;
using namespace std;

void testAccess() {
    cout<<"constructor"<<endl;
    MatrixDynamicSize<double> mat;
    mat.resize(2,10);

    for(uint32_t i=0;i<mat.size();i++) {
        mat(i)=i+1;
    }

    cout<<"size of mat=["<<mat.rows()<<" , "<<mat.cols()<<"]\n";

    Eigen::Map<Eigen::ArrayXXd> map(mat.data(),mat.rows(),mat.cols());
    //map.resize();

    cout<<"size of map=["<<map.rows()<<" , "<<map.cols()<<"]\n";

    cout<<map<<endl;
    cout<<"finished"<<endl;
}

void testCopy() {
    MatrixDynamicSize<size_t> matA,matB;

    matA.resize(3,4);
    size_t temp=0;
    for(auto & i : matA) {
        i=(temp++);
    }

    Eigen::Map<Eigen::Array<size_t,-1,-1>> map(matA.data(),matA.rows(),matA.cols());

    cout<<map<<endl;

    //auto funPtr=&std::array<void *,10>::size;

    matB=matA;

    cout<<Eigen::Map<Eigen::Array<size_t,-1,-1>>(matB.data(),matB.rows(),matB.cols())<<endl;

}

class Double
{
public:
    double val;
};

void testCustomTypes() {
    MatrixDynamicSize<MatrixFixedSize<Double,3,4>> matMat;
    matMat.resize(2,3);

}

void testLoop(uint32_t loopN) {
    MatrixDynamicSize<Double> mat;
    std::array<uint32_t,2> sizeList[]={{2,6},{10,2},{5,36},{6,30},{40,1},{1,100}};
    for(uint32_t i=0;i<loopN;i++) {
        uint32_t idx=i%(sizeof(sizeList)/sizeof(sizeList[0]));
        //cout<<"Loop "<<i<<", size=["<<sizeList[idx][0]<<" , "<<sizeList[idx][1]<<']'<<endl;
        mat.resize(sizeList[idx][0]*10,sizeList[idx][1]*10);
    }
}
