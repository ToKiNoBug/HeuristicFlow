#include <iostream>
#include "includes.h"
#include "lab4NSGA2.h"
#include "lab4NSGA3.h"

OptimT_MAKE_GLOBAL
using namespace std;
using namespace Eigen;
int main() {

    double c[]={M_PI_2,
                M_PI/3,
                M_PI/6,
                M_2_PI,
                M_PI*3,
                M_PI*4,
                M_PI*6,
                M_SQRT1_2,
                M_SQRT2
               };
    /*
1744352159520604142
458953659822443224
114738414981121388
993512999210214248
2111843339415065010
4223686678804044427
2111843339388979417
566232695896337324
1464156825348615605
*/

    for(int i=0;i<sizeof(c)/sizeof(c[0]);i++) {
        cout<<findCodeOf(c[i])<<endl;
    }

    return 0;

    testNSGA3Expri();


    return 0;
}
