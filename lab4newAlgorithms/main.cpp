#include <iostream>
#include "lab4NSGA2.h"
#include "lab4NSGA3.h"

OptimT_MAKE_GLOBAL
using namespace std;
int main() {

    Eigen::Matrix3d P;
    P<<1,0,0,
       0,2,0,
       0,0,3;
    {
        Eigen::Matrix3d temp;
        temp.leftCols(2)=P.rightCols(2);
        temp.col(2)=P.col(0);
        double r1=OptimT::randD(1e-9,10),r2=OptimT::randD(1e-9,10);
        auto temptemp=((temp*r1+P*r2)/(r1+r2)).eval();
        P=temptemp;

    }
    {
        Eigen::Matrix3d temp;
        temp.leftCols(2)=P.rightCols(2);
        temp.col(2)=P.col(0);
        double r1=OptimT::randD(1e-9,10),r2=OptimT::randD(1e-9,10);
        auto temptemp=((temp*r1+P*r2)/(r1+r2)).eval();
        P=temptemp;

    }

    Eigen::Array3d intercept=sample2Intercept(P);
    cout<<"Matrix P=\n"<<P<<endl<<endl;
    cout<<"intercept=\n"<<intercept<<endl;

    return 0;
}
