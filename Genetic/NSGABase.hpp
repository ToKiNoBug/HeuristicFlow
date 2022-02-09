/*
 Copyright © 2022  TokiNoBug
This file is part of OptimTemplates.

    OptimTemplates is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OptimTemplates is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OptimTemplates.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef OptimT_NSGABASE_HPP
#define OptimT_NSGABASE_HPP

#include "MOGABase.hpp"

#ifdef OptimT_NSGA_USE_THREADS
#ifndef OptimT_USE_THREADS
#error You allowed parallelize in NSGA2 but not in global.  \
    Macro OptimT_NSGA2_USE_THREADS can only be defined when OptimT_USE_THREADS is defined.
#endif
#endif

namespace OptimT {

template<typename Var_t,
        size_t ObjNum,
        typename Fitness_t,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class Args_t=void>
class NSGABase
    : public MOGABase<Var_t,ObjNum,Fitness_t,fOpt,rOpt,pfOpt,Args_t>
{
public:
    NSGABase() {};
    virtual ~NSGABase() {};

    using Base_t = MOGABase<Var_t,ObjNum,Fitness_t,fOpt,rOpt,pfOpt,Args_t>;
    OptimT_MAKE_GABASE_TYPES
    
    /**
     * @brief basical unit for NS
     * 
     */
    struct infoUnitBase
    {
    public:
        /** @brief Genes in population that strong domain this gene
        */
        size_t domainedByNum;

        /** @brief iterator to related gene
        */
        GeneIt_t iterator;
    };  //  infoUnitBase

protected:
    //calculate domainedByNum
    virtual void calculateDominatedNum(infoUnitBase ** pop,
        const size_t popSizeBefore) const {
#ifdef OptimT_NSGA_USE_THREADS
        static const size_t thN=OtGlobal::threadNum();
#pragma omp parallel for
        for(size_t begIdx=0;begIdx<thN;begIdx++) {

            for(size_t ed=begIdx;ed<popSizeBefore;ed+=thN) {
                pop[ed]->domainedByNum=0;
                for(size_t er=0;er<popSizeBefore;er++) {
                    if(er==ed)
                        continue;
                    pop[ed]->domainedByNum+=
                            Base_t::isStrongDomain(&(pop[er]->iterator->_Fitness),
                                           &(pop[ed]->iterator->_Fitness));
                }
            }
        }

#else
        for(size_t ed=0;ed<popSizeBefore;ed++) {
            pop[ed]->domainedByNum=0;
            for(size_t er=0;er<popSizeBefore;er++) {
                if(er==ed)
                    continue;
                pop[ed]->domainedByNum+=
                        Base_t::isStrongDomain(&(pop[er]->iterator->_Fitness),
                                       &(pop[ed]->iterator->_Fitness));
            }
        }
#endif

    }   //calculateDominatedNum()

    void updatePF(const infoUnitBase ** pfs,const size_t curFrontSize) {
            this->_pfGenes.clear();
            for(size_t i=0;i<curFrontSize;i++) {
                this->_pfGenes.emplace(&*(pfs[i]->iterator));
            }
            if(this->prevFrontSize!=curFrontSize) {
                this->_failTimes=0;
                this->prevFrontSize=curFrontSize;
            }
            else {
                size_t checkSum=this->makePFCheckSum();

                if(this->prevPFCheckSum==checkSum) {
                    this->_failTimes++;
                } else {
                    this->_failTimes=0;
                    this->prevPFCheckSum=checkSum;
                }
            }
    }   //updatePF()


};  //NSGABase


#define OptimT_MAKE_NSGABASE_TYPES \
OptimT_MAKE_GABASE_TYPES \
using infoUnitBase_t = typename Base_t::infoUnitBase;

}   // OptimT

#endif  //  NSGABASE_HPP
