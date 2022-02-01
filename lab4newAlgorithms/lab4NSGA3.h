#ifndef OptimT_LAB4NSGA3_H
#define OptimT_LAB4NSGA3_H

#include <Eigen/Dense>
#include <vector>

#include <OptimTemplates/Global>

std::vector<Eigen::ArrayXd> makeReferencePoints(const uint64_t dimN,const uint64_t precision);

#endif  //  OptimT_LAB4NSGA3_H