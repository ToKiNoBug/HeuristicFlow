#include "def_TestingPSO.h"
using namespace OptimT;
using namespace std;
void testPSOBase() {
    PSOBase<std::array<double,10>,
            double,
            RecordOption::DONT_RECORD_FITNESS> * noRec;

    PSOBase<std::array<double,10>,
            double,
            RecordOption::RECORD_FITNESS> * doRec;

    noRec=doRec;

    std::cout<<"sizeof noRec="<<sizeof(typeof (*noRec))<<std::endl;
    std::cout<<"sizeof doRec="<<sizeof(typeof (*doRec))<<std::endl;

    std::array<double,10> a;
    a[0]=0;

    PSO<std::array<double,10>,
            10,
            StdArray,
            FITNESS_GREATER_BETTER,
            DONT_RECORD_FITNESS> psoNoRec;

    PSO<std::array<double,10>,
            10,
            StdArray,
            FITNESS_GREATER_BETTER,
            RECORD_FITNESS> psoDoRec;

    noRec=&psoDoRec;
    doRec=&psoDoRec;

    noRec=&psoNoRec;

    ///doRec=&psoNoRec;
    /// code above failed to compile,
    /// because PSO without fitness recording is not derived from PSOBase with recording
    ///error: incompatible pointer types
    ///assigning to
    /// 'PSOBase<std::array<double, 10>, double, RecordOption::RECORD_FITNESS> *'
    /// from
    /// 'PSO<std::array<double, 10>, 10, StdArray, FITNESS_GREATER_BETTER, DONT_RECORD_FITNESS> *'

}
