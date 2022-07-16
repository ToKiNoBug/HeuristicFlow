#ifndef HEU_SOGASELECTORS_HPP
#define HEU_SOGASELECTORS_HPP

#include "GABase.hpp"

#include <unordered_map>

namespace heu {
namespace internal {
template <SelectMethod sm>
class SOGASelector {
  // static_assert(false);
};

template <>
class SOGASelector<SelectMethod::Truncation> {
 protected:
  template <class this_t>
  void __impl___impl_select() noexcept {
    HEU_MAKE_GABASE_TYPES(this_t)
    std::list<Gene_t>& popRef = static_cast<this_t*>(this)->_population;
    std::vector<GeneIt_t> sortSpace;
    sortSpace.clear();
    sortSpace.reserve((popRef.size()));

    for (GeneIt_t it = popRef.begin(); it != popRef.end(); ++it) {
      sortSpace.emplace_back(it);
    }

    std::sort(sortSpace.begin(), sortSpace.end(), this_t::GeneItCompareFun);

    while (static_cast<this_t*>(this)->_population.size() >
           static_cast<this_t*>(this)->_option.populationSize) {
      static_cast<this_t*>(this)->_population.erase(sortSpace.back());
      sortSpace.pop_back();
    }

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(sortSpace.front());
  }

 private:
};

template <>
class SOGASelector<SelectMethod::RouletteWheel> {
 protected:
  template <class this_t>
  void __impl___impl_select() noexcept {
    HEU_MAKE_GABASE_TYPES(this_t)

    /*
     * In this function, the linked list `eliminatedGenes` stores all candidates to be eliminated.
     * At first it has all candidates, and as the loop running, we erase a candidate from
     * `eliminatedGenes` where stores candidates waiting to be eliminated, symbolizing this
     * candidate is selected.
     *
     * When the size of `eliminatedGenes` equals to eliminate num, all candidates inside
     * this hashmap will be erased from the linked list `_population`
     */

    std::list<std::pair<GeneIt_t, double> >
        eliminatedGenes;  // Use the iterator as first and processed fitness as second

    // number of candidates that need to be eliminated.
    const int eliminateNum = static_cast<this_t*>(this)->_population.size() -
                             static_cast<this_t*>(this)->_option.populationSize;

    // A cache to store the best fitness value in the previous generator cause the best gene may be
    // eliminated
    const double previousBestFitness = static_cast<this_t*>(this)->bestFitness();

    // eliminatedGenes.reserve(static_cast<this_t*>(this)->_population.size());

    // fill eliminatedGenes
    for (GeneIt_t it = static_cast<this_t*>(this)->_population.begin();
         it != static_cast<this_t*>(this)->_population.end(); ++it) {
      // If fitness option is FITNESS_LESS_BETTER, take the inverse value
      if constexpr (this_t::FitnessOpt == FitnessOption::FITNESS_GREATER_BETTER) {
        eliminatedGenes.emplace_back(std::make_pair(it, it->_Fitness));
      } else {
        eliminatedGenes.emplace_back(std::make_pair(it, -it->_Fitness));
      }
    }

    // erase all selected candidates from eliminatedGenes
    while (eliminatedGenes.size() > eliminateNum) {
      double minFitness = eliminatedGenes.front().second;
      double processedFitnessSum = 0;
      // find the worst fitness
      for (auto& pair : eliminatedGenes) {
        minFitness = std::min(minFitness, pair.second);
      }

      // update and compute the sum of processed fitness
      for (auto& pair : eliminatedGenes) {
        pair.second -= minFitness;
        processedFitnessSum += pair.second;
      }

      // erase a solution from eliminatedGenes by its fitness
      if (processedFitnessSum > 1e-10) {
        double r = randD(0, processedFitnessSum);
        for (auto it = eliminatedGenes.begin(); it != eliminatedGenes.end(); ++it) {
          r -= it->second;
          if (r <= 0) {
            eliminatedGenes.erase(it);
            break;
          }
        }

        if (r > 0) {  // unlikely
          eliminatedGenes.pop_back();
        }

      } else {
        // select stochastically
        const int sizeBeforeErasement = eliminatedGenes.size();
        for (auto it = eliminatedGenes.begin(); it != eliminatedGenes.end(); ++it) {
          if (randIdx(sizeBeforeErasement) == 0) {
            eliminatedGenes.erase(it);
            break;
          }
        }
      }
    }

    // erase all eliminated candidates from the linked list
    for (auto& pair : eliminatedGenes) {
      static_cast<this_t*>(this)->_population.erase(pair.first);
    }

    GeneIt_t bestGeneIt = static_cast<this_t*>(this)->_population.begin();
    for (GeneIt_t it = static_cast<this_t*>(this)->_population.begin();
         it != static_cast<this_t*>(this)->_population.end(); ++it) {
      if (this_t::isBetter(it->_Fitness, bestGeneIt->_Fitness)) {
        bestGeneIt = it;
      }
    }

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(bestGeneIt, previousBestFitness);
  }
};

}  // namespace internal
}  // namespace heu

#endif  // HEU_SOGASELECTORS_HPP