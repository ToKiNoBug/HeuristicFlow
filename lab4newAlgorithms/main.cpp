#include <iostream>
#include "testNsga2.h"
#include "lab4NSGA3.h"

OptimT_MAKE_GLOBAL
using namespace std;
int main() {
    cout<<OptimT::fractorial<uint64_t>(10)<<endl;
    cout<<OptimT::NchooseK(6,3)<<endl;
    makeReferencePoints(3,10);
    return 0;
}
