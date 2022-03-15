#include "def_TestingBoxes.h"
#include <iostream>
using namespace Heu;
using namespace std;
void test_Box_double() {
    //50 dim square box in [0,1]
    Heu::template BoxNdS<50,DoubleVectorOption::Eigen,
            10,encode<0,1>::code,encode<1,1>::code> box;

    static constexpr size_t D=box.dimensions();
    cout<<box.Flag;
    static constexpr double minV=box.min();
    static constexpr double maxV=box.max();

    static const size_t szBox=sizeof(box);


}
