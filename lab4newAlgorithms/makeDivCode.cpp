#include "includes.h"
#include <cmath>
#include <iostream>

using namespace std;

DivCode findCodeOf(const double d,double * error) {
    if(d<=0||d==1)
        exit(114514);

    const bool greaterThan1=d>1;
    const int32_t beg=100000,end=0x7FFFFFFF-10;

    std::array<std::pair<DivCode,double>,4> result;

#pragma omp parallel for
    for(int32_t begVal=beg;begVal<beg+4;begVal++) {
        const int32_t id=begVal-beg;

        int32_t bestNum;
        uint32_t bestDen;
        double minError=1e10;
        if(greaterThan1)
            for(int32_t num=begVal;num<end;num+=4) {
                uint32_t den=std::round(num/d);
                double error=std::abs(double(num)/double(den)-d);
                if(minError>error) {
                    minError=error;
                    bestNum=num;
                    bestDen=den;
                }
            }
        else
            for(uint32_t den=begVal;den<end;den+=4) {
                int32_t num=std::round(den*d);
                double error=abs(double(num)/double(den)-d);
                if(minError>error) {
                    minError=error;
                    bestNum=num;
                    bestDen=den;
                }
            }


        result[id].first=encodeFun(bestNum,bestDen);
        result[id].second=minError;

    }
    int minIdx=0;
    for(int i=0;i<4;i++) {
        //cout<<"Code "<<i<<" = "<<result[i].first<<" , error = "<<log10(result[i].second+1e-300)<<endl;
        if(result[i].second<result[minIdx].second) {
            minIdx=i;
        }
    }

    if(error!=nullptr) {
        *error=result[minIdx].second;
    }
    return result[minIdx].first;
}
