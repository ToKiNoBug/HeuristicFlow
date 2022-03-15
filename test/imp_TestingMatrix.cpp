/*
 Copyright © 2022  TokiNoBug
This file is part of HeuristicFlow.

    HeuristicFlow is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HeuristicFlow is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HeuristicFlow.  If not, see <https://www.gnu.org/licenses/>.

*/

#include "def_TestingMatrix.h"
#include <HeuristicFlow/SimpleMatrix>
#include <stdint.h>
#include <iostream>
#include <Eigen/Dense>
#include <array>
#include <ctime>
using namespace Heu;
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
    clock_t c=clock();
    for(uint32_t i=0;i<loopN;i++) {
        uint32_t idx=i%(sizeof(sizeList)/sizeof(sizeList[0]));
        //cout<<"Loop "<<i<<", size=["<<sizeList[idx][0]<<" , "<<sizeList[idx][1]<<']'<<endl;
        mat.resize(sizeList[idx][0]*10,sizeList[idx][1]*10);
    }
    cout<<double(clock()-c)*1000/CLOCKS_PER_SEC/loopN*10e3<<" ms per allocate"<<endl;
}


void testInverse() {
    static const size_t N=Runtime;
    MatrixDynamicSize<double> A(50,50),iA(50,50);
    
    for(auto & i : A) {
        i=randD(-1,1);
    }

    cout<<"A=[";
    for(size_t c=0;c<A.cols();c++) {
        for(size_t r=0;r<A.rows();r++) {
            cout<<A(r,c)<<',';
        }
        cout<<";\n";
    }
    cout<<"];\n\n\n\n";

    InverseMatrix_LU<double,N>(A,&iA);

    cout<<"iA=[";
    for(size_t c=0;c<iA.cols();c++) {
        for(size_t r=0;r<iA.rows();r++) {
            cout<<iA(r,c)<<',';
        }
        cout<<";\n";
    }
    cout<<"];\n\n\n\n"<<endl;

}

void testProduct() {
    MatrixDynamicSize<double> A(10,6);
    MatrixDynamicSize<double> B(6,4);
    MatrixDynamicSize<double> C;
    for(auto & i : A) {
        i=randD(-2,2);
    }
    for(auto & i : B) {
        i=randD(-1,1);
    }

    {
    Eigen::Map<Eigen::ArrayXXd> mA(A.data(),A.rows(),A.cols()),
            mB(B.data(),B.rows(),B.cols());

    cout<<"A=["<<mA<<"];\n\n\n"<<endl;
    cout<<"B=["<<mB<<"];\n\n\n"<<endl;
    }

    MatrixProduct(A,B,&C);

    cout<<"C=["<<Eigen::Map<Eigen::ArrayXXd>(C.data(),C.rows(),C.cols())<<"];\n\n\n"<<endl;
}
