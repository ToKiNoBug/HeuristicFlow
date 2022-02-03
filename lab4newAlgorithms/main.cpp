#include <iostream>
#include "lab4NSGA2.h"
#include "lab4NSGA3.h"

OptimT_MAKE_GLOBAL
using namespace std;
using namespace Eigen;
int main() {
    /*
    Matrix3d P;
    P<<1,0,0,
       0,2,0,
       0,0,3;
    {
        Matrix3d temp;
        temp.leftCols(2)=P.rightCols(2);
        temp.col(2)=P.col(0);
        double r1=OptimT::randD(1e-9,10),r2=OptimT::randD(1e-9,10);
        auto temptemp=((temp*r1+P*r2)/(r1+r2)).eval();
        P=temptemp;

    }
    {
        Matrix3d temp;
        temp.leftCols(2)=P.rightCols(2);
        temp.col(2)=P.col(0);
        double r1=OptimT::randD(1e-9,10),r2=OptimT::randD(1e-9,10);
        auto temptemp=((temp*r1+P*r2)/(r1+r2)).eval();
        P=temptemp;

    }

    Array3d intercept=sample2Intercept(P);
    cout<<"Matrix P=\n"<<P<<endl<<endl;
    cout<<"intercept=\n"<<intercept<<endl;
    */

    Array<double,3,60> w;
    Array3d s;
    w.setRandom();
    s.setRandom();
    cout<<"w=[\n"<<w<<"];"<<endl;
    cout<<"s=["<<s.transpose()<<"]';"<<endl;

    auto wT_s=w.matrix().transpose()*s.matrix();
    auto wT_s_w=w.rowwise()*(wT_s.array().transpose());
    typeof(w) norm_wTsw=wT_s_w.rowwise()/(w.colwise().squaredNorm());
    auto s_sub_norm_wTsw=norm_wTsw.colwise()-s;
    auto distance=s_sub_norm_wTsw.colwise().squaredNorm();
    cout<<"distance=["<<distance<<"];"<<endl;




    return 0;
}
