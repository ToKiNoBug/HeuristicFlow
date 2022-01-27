#ifndef MOGAABSTRACT_HPP
#define MOGAABSTRACT_HPP

#include "./GABase.hpp"
#include <queue>
#include <unordered_set>

namespace OptimT {
    
///whether to protect pareto front when mutation or not
enum PFOption : unsigned char {
    PARETO_FRONT_DONT_MUTATE=true,
    PARETO_FRONT_CAN_MUTATE=false
};

/**
   *  @brief Base class for multi-objective genetic algorithm solver.
   *
   *  @tparam Var_t  Type of decisition variable.
   *  @tparam ObjNum Numbers of objectives.
   *  @tparam Fitness_t Type of fitness value.
   *  @tparam fOpt Whether greater fitness value means better.
   *  @tparam rOpt Whether the solver records fitness changelog.
   *  @tparam pfOpt Whether to protect the Pareto front from mutation.
   *  @tparam ...Args Type of other parameters.
  */
template<typename Var_t,
        size_t ObjNum,
        typename Fitness_t,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class ...Args>
class MOGAAbstract 
    : public GABase<Var_t,Fitness_t,rOpt,Args...>
{
public:
    MOGAAbstract()  {};
    virtual ~MOGAAbstract() {};

    using Base_t = GABase<Var_t,Fitness_t,rOpt,Args...>;
    OptimT_MAKE_GABASE_TYPES

    ///get pareto front in vec
    void paretoFront(std::vector<Fitness_t> & front) const {
        front.clear();  front.reserve(_pfGenes.size());
        for(const Gene* i : _pfGenes) {
            front.emplace_back(i->_Fitness);
        }
        return;
    }

    void paretoFront(std::vector<std::pair<const Var_t*,const Fitness_t*>> & front) const {
        front.clear();
        front.reserve(_pfGenes.size());
        for(const Gene* i : _pfGenes) {
            front.emplace_back(std::make_pair(&(i->self),&(i->_Fitness)));
        }
    }

    const std::unordered_set<const Gene*> & pfGenes() const {
        return _pfGenes;
    }


protected:
    size_t prevFrontSize;
    size_t prevPFCheckSum;
    std::unordered_set<const Gene*> _pfGenes;

    virtual size_t makePFCheckSum() const {
        std::vector<const Gene*> pfvec;
        pfvec.reserve(_pfGenes.size());
        for(auto i : _pfGenes) {
            pfvec.emplace_back(i);
        }
        std::sort(pfvec.begin(),pfvec.end());

        static const auto hashFun=_pfGenes.hash_function();
        size_t checkSum=hashFun(pfvec.front());
        for(size_t i=1;i<pfvec.size();i++) {
            checkSum^=hashFun(pfvec[i]);
        }
        return checkSum;
    }

    ///mutate operation
    virtual void mutate() {
        for(auto it=Base_t::_population.begin();it!=Base_t::_population.end();++it) {
            if(randD()<=Base_t::_option.mutateProb) {
                if(pfOpt){
                    if(_pfGenes.find(&*it)!=_pfGenes.end()) {
                        continue;
                    }
                }
                Base_t::_mutateFun(&it->self,&Base_t::args());
                it->setUncalculated();
            }
        }
    }

private:
    
#ifndef OptimT_NO_STATICASSERT
    static_assert(std::integral_constant<bool,(ObjNum!=1)>::value,
    "OptimTemplates : You used less than 1 objective in NSGA2");
#endif

};

}


#endif //   MOGAABSTRACT_HPP
