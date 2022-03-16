#include "def_TestingBoxes.h"
#include <iostream>
using namespace Heu;
using namespace std;

void test_Box_double() {
    BoxNdS<50,DoubleVectorOption::Eigen,
            10,encode<0,1>::code,encode<1,1>::code,encode<1,50>::code> box0;
    //50 dim square box in [0,1]
    BoxXdN<DoubleVectorOption::Eigen> box;

    //cout<<box.Flag<<endl;

    box.max().setConstant(50,1,1.0);
    box.min().setConstant(50,1,0.0);

    box0.min();

    cout<<box0.dimensions()<<endl;
    cout<<box0.learnRate()<<endl;
    cout<<"sizeof(box0) = "<<sizeof(box0)<<endl;

}

void test_Box_bool() {
    BoxNb<10> BNb;
    BoxXb<> BXb;

    BXb.setDimensions(400);

    BNb.max();
    BNb.min();

    BXb.max();
    BXb.min();

    cout<<BNb.dimensions();
    cout<<BXb.dimensions();

    cout<<"sizeof BNb="<<sizeof(BNb)<<endl;
    cout<<"sizeof BXb="<<sizeof(BXb)<<endl;

}
