#ifndef OptimT_LAB4NSGA3_H
#define OptimT_LAB4NSGA3_H

#include "includes.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>

Eigen::ArrayXd sample2Intercept(Eigen::MatrixXd);
std::vector<Eigen::ArrayXd> makeReferencePoints(const uint64_t dimN,const uint64_t precision);

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
        _precision=4;
    }

    virtual ~testNSGA3() {

    }

    using Base_t = OptimT::NSGABase<Eigen::Array2d,
        3,
        Eigen::Array3d,
        OptimT::FITNESS_LESS_BETTER,
        OptimT::RECORD_FITNESS,
        OptimT::PARETO_FRONT_DONT_MUTATE>;

    OptimT_MAKE_NSGABASE_TYPES
    using RefPoint = size_t;
    /*
    struct RefPoint
    {
    public:
        size_t nicheCount;
    };
    */

    struct infoUnit3 : public infoUnitBase_t
    {
    public:
        Eigen::Array3d translatedFitness;
        size_t closestIdx;
        double distance;
    };

    inline size_t precision() const {
        return _precision;
    }

    inline void setPrecision(size_t p) {
        _precision=p;
    }

protected:
    size_t _precision;
    Eigen::Array<double,3,Eigen::Dynamic> referencePoses;

    void customOptWhenInitialization() {
        auto x=makeReferencePoints(ObjNum,_precision);
        referencePoses.resize(3,x.size());
        for(size_t i=0;i<x.size();i++) {
            referencePoses.col(i)=x[i];
        }
    }

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
            pop.back().closestIdx=-1;
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

        std::unordered_set<infoUnit3*> selected;
        selected.reserve(this->_option.populationSize);

        std::vector<infoUnit3*> * FlPtr=nullptr;
        bool needRefPoint=false;

        while(true) {
            if(selected.size()==this->_option.populationSize) {
                needRefPoint=false;
                break;
            }
            if(selected.size()+pfLayers.front().size()>this->_option.populationSize) {
                needRefPoint=true;
                FlPtr=&pfLayers.front();
                break;
            }

            for(infoUnit3* i : pfLayers.front()) {
                selected.emplace(i);
            }
            pfLayers.pop_front();

        }

        if(needRefPoint) {
            ///Normalize procedure
            std::unordered_set<infoUnit3*> Fl;
            Fl.reserve(FlPtr->size());
            for(auto i : *FlPtr) {
                Fl.emplace(i);
            }
            normalize(selected,Fl);
            std::unordered_map<size_t,RefPoint> refPoints;
            refPoints.reserve(referencePoses.cols());
            for(size_t i=0;i<referencePoses.cols();i++) {
                refPoints[i]=0;
            }
            ///Associate procedure
            associate(selected);
            associate(Fl);
            nichePreservation(&selected,&Fl,&refPoints);
        }

        //erase all unselected genes
        for(auto i : sortSpace) {
            if(selected.find(i)==selected.end()) {
                this->_population.erase(i->iterator);
            }
        }

    }   //  end select

    void normalize(const std::unordered_set<infoUnit3*> & selected,
        const std::unordered_set<infoUnit3*> & Fl) const {
        auto it=(selected.begin());
        Eigen::Array3d ideal=(*it)->iterator->_Fitness;
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

        Eigen::Array3d intercept;
        
        extremePoints2Intercept(extremePoints,intercept);

        for(auto i : selected) {
            i->translatedFitness=(i->iterator->_Fitness-ideal)/intercept;
        }
        for(auto i : Fl) {
            i->translatedFitness=(i->iterator->_Fitness-ideal)/intercept;
        }
        
    }

    void associate(const std::unordered_set<infoUnit3*> & st) const {

        for(auto i : st) {
            const auto & w=referencePoses;
            const auto & s=i->translatedFitness;

            auto wT_s=w.matrix().transpose()*s.matrix();
            auto wT_s_w=w.rowwise()*(wT_s.array().transpose());
            Eigen::Array3Xd norm_wTsw=wT_s_w.rowwise()/(w.colwise().squaredNorm());
            auto s_sub_norm_wTsw=norm_wTsw.colwise()-s;
            auto distance=s_sub_norm_wTsw.colwise().squaredNorm();

            int minDistanceIdx;
            i->distance=distance.minCoeff(&minDistanceIdx);
            i->closestIdx=minDistanceIdx;
        }
    }

    void nichePreservation(std::unordered_set<infoUnit3*> * selected,
            std::unordered_set<infoUnit3*> * Fl,
            std::unordered_map<size_t,RefPoint> * refPoints) {
        for(auto i : *selected) {
            refPoints->operator[](i->closestIdx)++;
        }
        std::vector<std::unordered_map<size_t,RefPoint>::iterator> minNicheIterators;
        minNicheIterators.reserve(refPoints->size());

        std::vector<infoUnit3*> associatedGenesInFl;
        associatedGenesInFl.reserve(Fl->size());

        while(selected->size()<=this->_option.populationSize) {
            findMinSet(*refPoints,&minNicheIterators);
            auto curRefPoint=minNicheIterators[size_t(OptimT::randD(0,minNicheIterators.size()))];
            size_t rhoJ=curRefPoint->second;
            findAssociated(*Fl,curRefPoint,&associatedGenesInFl);


            if(!associatedGenesInFl.empty()) {
                infoUnit3 * pickedGene=nullptr;
                if(rhoJ==0) {
                    //find element in associatedGenesInFl with minimum distance
                    infoUnit3 * minGene=associatedGenesInFl.front();
                    for(auto i : associatedGenesInFl) {
                        if(i->distance<minGene->distance) {
                            minGene=i;
                        }
                    }
                    pickedGene=minGene;
                }
                else {
                    //pick a random member in associatedGenesInFl
                    pickedGene=associatedGenesInFl[size_t(OptimT::randD(0,associatedGenesInFl.size()))];
                }
                selected->emplace(pickedGene);
                Fl->erase(pickedGene);
                curRefPoint->second++;
            }
            else {
                refPoints->erase(curRefPoint);
            }
        }   //  end while
    }

    inline static void findMinSet(std::unordered_map<size_t,RefPoint> & refPoints,
            std::vector<std::unordered_map<size_t,RefPoint>::iterator> * minNicheIterators) {
        minNicheIterators->clear();
        size_t minNiche=-1;
        for(auto i : refPoints) {
            minNiche=std::min(minNiche,i.second);
        }
        for(auto it=refPoints.begin();it!=refPoints.end();++it) {
            if(it->second==minNiche) {
                minNicheIterators->emplace_back(it);
            }
        }
    }

    inline static void findAssociated(const std::unordered_set<infoUnit3*> & Fl,
            const std::unordered_map<size_t,RefPoint>::iterator & refP,
            std::vector<infoUnit3*> * associatedGenesInFl) {
        associatedGenesInFl->clear();

        for(auto i : Fl) {
            if(i->closestIdx==refP->first) {
                associatedGenesInFl->emplace_back(i);
            }
        }

    }

    inline static void extremePoints2Intercept(const Eigen::Array33d & P,Eigen::Array3d & intercept) {
        auto P_transpose_inv=P.transpose().matrix().inverse();
        auto ONE=Eigen::Matrix<double,Eigen::Dynamic,1>::Ones(P.cols(),1);
        auto one_div_intercept=(P_transpose_inv*ONE).array();
        intercept=1.0/one_div_intercept;
    }
};


#endif  //  OptimT_LAB4NSGA3_H
