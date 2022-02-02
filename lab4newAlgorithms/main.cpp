#include <iostream>
#include "lab4NSGA2.h"
#include "lab4NSGA3.h"

OptimT_MAKE_GLOBAL
using namespace std;
int main() {

    Eigen::Matrix3d P;
    P.col(0)=Eigen::Array3d({1,0,0});
    P.col(1)=Eigen::Array3d({0,2,0});
    P.col(2)=Eigen::Array3d({0,0,3});

    Eigen::Array3d intercept=sample2Intercept(P);

    cout<<"intercept=\n"<<intercept<<endl;

    return 0;
}
