#ifndef testNsga3_H
#define testNsga3_H

#include <Eigen/Dense>

#define OptimT_NO_OUTPUT
#define OptimT_DO_PARRALLELIZE
#include <OptimTemplates/Global>
#include <OptimTemplates/Genetic>

const size_t spaceDim=3;
const size_t targetNum=6;

class testNsga3 
    : public OptimT::GABase<
        Eigen::Array<double,spaceDim,1>,
        Eigen::Array<double,targetNum,1>,
        OptimT::RecordOption::RECORD_FITNESS>
{
public:
    testNsga3() {};
    ~testNsga3() {};
    using Base_t = OptimT::GABase<
        Eigen::Array<double,spaceDim,1>,
        Eigen::Array<double,targetNum,1>,
        OptimT::RecordOption::RECORD_FITNESS>;
        
    OPTIMT_MAKE_GABASE_TYPES

};

#endif