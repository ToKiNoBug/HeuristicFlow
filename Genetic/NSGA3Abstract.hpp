#ifndef OptimT_NSGA3ABSTRACT_HPP
#define OptimT_NSGA3ABSTRACT_HPP

#include "NSGABase.hpp"

namespace OptimT {


template<typename Var_t,
        size_t ObjNum,
        DoubleVectorOption DVO,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class Args_t>
class NSGA3Abstract
    : public NSGABase<Var_t,ObjNum,FitnessVec_t<DVO,ObjNum>,fOpt,rOpt,pfOpt,Args_t>
{
    NSGA3Base() {};
    virtual ~NSGA3Base() {};
    using Base_t = NSGABase<Var_t,ObjNum,FitnessVec_t<DVO,ObjNum>,fOpt,rOpt,pfOpt,Args_t>;
    OptimT_MAKE_NSGABASE_TYPES
    using Fitness_t = FitnessVec_t<DVO,ObjNum>;
    static_assert(DVO!=DoubleVectorOption::Custom,
        "Using custom double container as fitness isn't supported");
#ifdef EIGEN_CORE_H
    using RefMat_t= typename std::conditional<DVO==Std,
        MatrixDynamicSize<double>,
        typename std::conditional<
            ObjNum==Dynamic,
            Eigen::ArrayXXd,
            Eigen::<double,ObjNum,Eigen::Dynamic>
        ::type>::type;
#else
    using RefMat_t = MatrixDynamicSize<double>;
    static_assert(DVO!=DoubleVectorOption::Eigen,
        "Include Eigen before using Eigen arrays as Fitness types");
#endif

const RefMat_t & referencePoints() const {
    return referencePoses;
}

struct infoUnit3 : public infoUnitBase
{
    Fitness_t translatedFitness;
    size_t cloestRefPoint;
    double distance;
};

protected:
    RefMat_t referencePoses;

    void computeReferencePointPoses(const size_t dimN,
        const size_t precision,
        std::vector<Fitness_t> * dst) const {

        if(precision<=0) {
            exit(114514);
        }

        std::vector<Fitness_t> points;
        
        points.reserve(OptimT::NchooseK(dimN+precision-1,precision));

        pri_startRP(dimN,precision,&points);
    }

private:
    void pri_makeRP(const size_t dimN,const size_t precision,
        const size_t curDim,const size_t curP,const size_t accum,
        Fitness_t * rec,
        std::vector<Fitness_t> * dst) const {

        if(curDim+1>=dimN) {
            rec->operator[](dimN-1)=1.0-double(accum)/precision;
            dst->emplace_back(*rec);
            return;
        }
        

        for(size_t p=0;p+accum<=precision;p++) {
            if(curDim>=0)
                rec->operator[](curDim)=double(p)/precision;
            pri_makeRP(dimN,precision,curDim+1,p,accum+p,rec,dst);
        }
    }

    void pri_startRP(const size_t dimN,
        const size_t precision,
        vector<Eigen::ArrayXd> * dst) const {

        Fitness_t rec;

        if constexpr(ObjNum==Dynamic) {
            rec.resize(this->objectiveNum());
        }

        pri_makeRP(dimN,precision,0,0,0,&rec,dst);
    }
};

#define OptimT_MAKE_NSGA3ABSTRACT_TYPES \
OptimT_MAKE_NSGABASE_TYPES \
using RefMat_t = Base_t::RefMat_t; \
using Fitness_t = Base_t::Fitness_t; \
using infoUnit3 = Base_t::infoUnit3;

}

#endif  //  OptimT_NSGA3ABSTRACT_HPP