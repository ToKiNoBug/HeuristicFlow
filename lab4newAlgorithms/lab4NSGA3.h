#ifndef OptimT_LAB4NSGA3_H
#define OptimT_LAB4NSGA3_H

#include "includes.h"
#include <vector>

Eigen::ArrayXd sample2Intercept(Eigen::MatrixXd);

static const size_t VarDim=2;
static const size_t ObjNum=3;
/**
 * @brief Viennet function
 * 
 */
class testNSGA3
    : public OptimT::NSGABase<Eigen::Array2d,
        3,
        Eigen::Array3d,
        OptimT::FITNESS_LESS_BETTER,
        OptimT::RECORD_FITNESS,
        OptimT::PARETO_FRONT_DONT_MUTATE>
{
public:
    testNSGA3() {

    }

    ~testNSGA3() {

    }

    using Base_t = OptimT::NSGABase<Eigen::Array2d,
        3,
        Eigen::Array3d,
        OptimT::FITNESS_LESS_BETTER,
        OptimT::RECORD_FITNESS,
        OptimT::PARETO_FRONT_DONT_MUTATE>;

    OptimT_MAKE_NSGABASE_TYPES

    struct RefPoint
    {
    public:
        Eigen::Array2d Pos;
        size_t nicheCount;
    };

    struct infoUnit3 : public infoUnitBase_t
    {
    public:
        Eigen::Array3d translatedFitness;
        RefPoint * closest;
    };

protected:

    static bool sortByDominatedNum(const infoUnit3 * A,const infoUnit3* B) {
        if(A==B)
            return false;
        return (A->domainedByNum)<(B->domainedByNum);
    }

    void select() {
        const size_t popSizeBef=this->_population.size();
        std::vector<infoUnit3> pop;
        pop.reserve(popSizeBef);
        for(auto it=this->_population.begin();it!=this->_population.end();++it) {
            pop.emplace_back();
            pop.back().iterator=it;
            pop.back().closest=nullptr;
            pop.back().isSelected=false;
        }

        std::vector<infoUnit3*> sortSpace(popSizeBef);
        for(size_t i=0;i<popSizeBef;i++) {
            sortSpace[i]=pop.data()+i;
        }

        this->calculateDominatedNum((infoUnitBase_t**)sortSpace.data(),popSizeBef);

        std::list<std::vector<infoUnit3*>> pfLayers;
        std::sort(sortSpace.begin(),sortSpace.end(),sortByDominatedNum);

        size_t curDM=-1;
        for(auto i : sortSpace) {
            if(curDM!=i->domainedByNum) {
                curDM=i->domainedByNum;
                pfLayers.emplace_back();
                pfLayers.back().reserve(popSizeBef);
            }
            pfLayers.back().emplace_back(i);
        }

        this->updatePF((const infoUnitBase_t **)pfLayers.front().data(),pfLayers.front().size());

        std::vector<infoUnit3*> selected;
        selected.reserve(this->_option.populationSize);

        std::vector<infoUnit3*> * Fl=nullptr;
        bool needRefPoint=false;

        while(true) {
            if(selected.size()==this->_option.populationSize) {
                needRefPoint=false;
                break;
            }
            if(selected.size()+pfLayers.front().size()>this->_option.populationSize) {
                needRefPoint=true;
                Fl=&pfLayers.front();
                break;
            }

            for(infoUnit3* i : pfLayers.front()) {
                selected.emplace_back(i);
            }
            pfLayers.pop_front();

        }

        if(needRefPoint) {
            
        }



    }   //  end select

    void normalize(const std::vector<infoUnit3*> & selected,
        const std::vector<infoUnit3*> & Fl) const {
        Eigen::Array3d ideal=selected.front()->iterator->_Fitness;
        Eigen::Array33d extremePoints;
        extremePoints.colwise()=ideal;
        for(auto i : selected) {
            ideal=ideal.min(i->iterator->_Fitness);
            for(size_t objIdx=0;objIdx<ObjNum;objIdx++) {
                if(i->iterator->_Fitness[objIdx]>extremePoints(objIdx,objIdx)) {
                    extremePoints.col(objIdx)=i->iterator->_Fitness;
                }
            }
        }
        for(auto i : Fl) {
            ideal=ideal.min(i->iterator->_Fitness);
            for(size_t objIdx=0;objIdx<ObjNum;objIdx++) {
                if(i->iterator->_Fitness[objIdx]>extremePoints(objIdx,objIdx)) {
                    extremePoints.col(objIdx)=i->iterator->_Fitness;
                }
            }
        }

        extremePoints.colwise()-=ideal;

        for(auto i : selected) {
            i->translatedFitness=i->iterator->_Fitness-ideal;
        }
        for(auto i : Fl) {
            i->translatedFitness=i->iterator->_Fitness-ideal;
        }

        Eigen::Array3d intercept;

        auto P_transpose_inv=extremePoints.transpose().matrix().inverse();
        static const Eigen::Matrix<double,1,3> ONE13=Eigen::Matrix<double,1,3>::Ones();
        auto one_div_intercept=(ONE13*P_transpose_inv).array();

        intercept=1.0/one_div_intercept;

        
    }
};


std::vector<Eigen::ArrayXd> makeReferencePoints(const uint64_t dimN,const uint64_t precision);

#endif  //  OptimT_LAB4NSGA3_H