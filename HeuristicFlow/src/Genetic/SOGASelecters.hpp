#ifndef HEU_SOGASELECTORS_HPP
#define HEU_SOGASELECTORS_HPP

#include <unordered_map>

#include "GABase.hpp"

#include <iostream>

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

    std::list<std::pair<GeneIt_t, double>>
        eliminatedGenes;  // Use the iterator as first and processed fitness as second

    // number of candidates that need to be eliminated.
    const int eliminateNum = static_cast<this_t*>(this)->_population.size() -
                             static_cast<this_t*>(this)->_option.populationSize;

    // A cache to store the best fitness value in the previous generator cause the best gene may be
    // eliminated
    const double previousBestFitness = static_cast<this_t*>(this)->bestFitness();

    if (eliminateNum <= 0) {
      static_cast<this_t*>(this)->updateFailTimesAndBestGene(
          static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
      return;
    }
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
      if (processedFitnessSum > 0) {
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

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(
        static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
  }
};

template <>
class SOGASelector<SelectMethod::Tournament> {
 public:
  SOGASelector() { _tournamentSize = 3; }

  inline int& tournamentSize() noexcept { return _tournamentSize; }
  inline int tournamentSize() const noexcept { return _tournamentSize; };

  inline void setTournamentSize(int TS) noexcept {
    const bool tournament_size_should_be_greater_than_or_euqal_to_2 = TS >= 2;
    assert(tournament_size_should_be_greater_than_or_euqal_to_2);
    _tournamentSize = TS;
  }

 protected:
  int _tournamentSize;

 protected:
  template <class this_t>
  void __impl___impl_select() noexcept {
    HEU_MAKE_GABASE_TYPES(this_t);

    {
      const bool tournament_size_should_be_less_than_the_population_size =
          _tournamentSize < static_cast<this_t*>(this)->_option.populationSize;
      assert(tournament_size_should_be_less_than_the_population_size);
    }

    const int prevPopSize = static_cast<this_t*>(this)->_population.size();
    const double previousBestFitness = static_cast<this_t*>(this)->_bestGene->_Fitness;

    if (prevPopSize <= static_cast<this_t*>(this)->_option.populationSize) {
      static_cast<this_t*>(this)->updateFailTimesAndBestGene(
          static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
      return;
    }

    // the left _tournamentSize th elements are to be considered to be inside the tournament, while
    // the rest are outsize the tournament
    std::vector<Gene_t*> tournamentSpace(0);
    tournamentSpace.reserve(prevPopSize);

    std::unordered_map<Gene_t*, int, std::hash<void*>> selectCounter;
    selectCounter.reserve(prevPopSize);

    // fill tournamentSpace and selectCounter
    for (Gene_t& it : static_cast<this_t*>(this)->_population) {
      tournamentSpace.emplace_back(&it);
      selectCounter.emplace(&it, 0);
    }

    std::shuffle(tournamentSpace.begin(), tournamentSpace.end(), global_mt19937());

    // apply tournament selection
    for (int playTimes = 0; playTimes < static_cast<this_t*>(this)->_option.populationSize;
         playTimes++) {
      //   find the best gene inside the tournament
      Gene_t* bestGenePtr = tournamentSpace.front();
      for (int idx = 1; idx < _tournamentSize; idx++) {
        if (this_t::isBetter(tournamentSpace[idx]->_Fitness, bestGenePtr->_Fitness)) {
          bestGenePtr = tournamentSpace[idx];
        }
      }

      selectCounter[bestGenePtr]++;

      // pick three new candidates into the tournament
      for (int idxA = 0; idxA < _tournamentSize; idxA++) {
        const int idxB = randIdx(_tournamentSize, prevPopSize);
        std::swap<Gene_t*>(tournamentSpace[idxA], tournamentSpace[idxB]);
      }
    }

    // erase eliminated genes from population and selectCounter
    for (GeneIt_t it = static_cast<this_t*>(this)->_population.begin();
         it != static_cast<this_t*>(this)->_population.end();) {
      if (selectCounter[&*it] <= 0) {
        // This gene is not selected for even once. It is eliminated.
        selectCounter.erase(&*it);
        it = static_cast<this_t*>(this)->_population.erase(it);
        continue;
      }

      selectCounter[&*it]--;
      ++it;
    }

    // Now all eliminated genes are erased, and genes that are selected repeatedly hasn't been
    // copied. Searching at now will be slightly faster.

    GeneIt_t curBestGene = static_cast<this_t*>(this)->_population.begin();
    for (GeneIt_t it = static_cast<this_t*>(this)->_population.begin();
         it != static_cast<this_t*>(this)->_population.end(); ++it) {
      if (this_t::isBetter(it->_Fitness, curBestGene->_Fitness)) {
        curBestGene = it;
      }
    }

    for (auto& pair : selectCounter) {
      while (pair.second > 0) {
        // copy if a gene is selected for more than once
        static_cast<this_t*>(this)->_population.emplace_back(*pair.first);
        pair.second--;
      }
    }

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(curBestGene, previousBestFitness);

    // exit(114514);
  }
};

template <>
class SOGASelector<SelectMethod::MonteCarlo> {
 public:
  SOGASelector() = default;

 protected:
  template <class this_t>
  void __impl___impl_select() noexcept {
    HEU_MAKE_GABASE_TYPES(this_t);
    const double previousBestFitness = static_cast<this_t*>(this)->_bestGene->_Fitness;

    const int popSizeBeforeSelection = static_cast<this_t*>(this)->_population.size();

    if (popSizeBeforeSelection <= static_cast<this_t*>(this)->_option.populationSize) {
      static_cast<this_t*>(this)->updateFailTimesAndBestGene(
          static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
      return;
    }

    std::vector<GeneIt_t> iterators(0);
    iterators.reserve(popSizeBeforeSelection);

    for (GeneIt_t it = static_cast<this_t*>(this)->_population.begin();
         it != static_cast<this_t*>(this)->_population.end(); ++it) {
      iterators.emplace_back(it);
    }

    std::shuffle(iterators.begin(), iterators.end(), global_mt19937());

    const int eliminateNum =
        popSizeBeforeSelection - static_cast<this_t*>(this)->_option.populationSize;
    for (int idx = 0; idx < eliminateNum; idx++) {
      static_cast<this_t*>(this)->_population.erase(iterators[idx]);
    }

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(
        static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
  }
};

template <>
class SOGASelector<SelectMethod::Probability> {
 protected:
  template <class this_t>
  void __impl___impl_select() noexcept {
    HEU_MAKE_GABASE_TYPES(this_t)

    std::list<std::pair<GeneIt_t, double>>
        eliminatedGenes;  // Use the iterator as first and processed fitness as second
    const int prevPopSize = static_cast<this_t*>(this)->_population.size();
    const int K = static_cast<this_t*>(this)->_option.populationSize;
    // number of candidates that need to be eliminated.
    const int eliminateNum = prevPopSize - K;

    // A cache to store the best fitness value in the previous generator cause the best gene may be
    // eliminated
    const double previousBestFitness = static_cast<this_t*>(this)->bestFitness();

    if (eliminateNum <= 0) {
      static_cast<this_t*>(this)->updateFailTimesAndBestGene(
          static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
      return;
    }
    // assert(eliminateNum > 0);

    // fill fitness value and  Compute the probability in the RouttleWheel method

    {
      double minFitness = pinfD;
      for (GeneIt_t it = static_cast<this_t*>(this)->_population.begin();
           it != static_cast<this_t*>(this)->_population.end(); ++it) {
        // If fitness option is FITNESS_LESS_BETTER, take the inverse value
        if constexpr (this_t::FitnessOpt == FitnessOption::FITNESS_GREATER_BETTER) {
          eliminatedGenes.emplace_back(std::make_pair(it, it->_Fitness));
        } else {
          eliminatedGenes.emplace_back(std::make_pair(it, -it->_Fitness));
        }

        minFitness = std::min(eliminatedGenes.back().second, minFitness);
      }

      //
      double fitnessSum = 0;
      for (auto& pair : eliminatedGenes) {
        pair.second -= minFitness;
        fitnessSum += pair.second;
      }

      //
      if (fitnessSum <= 0) {
        for (auto& pair : eliminatedGenes) {
          pair.second = 1.0 / prevPopSize;
        }
      } else {
        for (auto& pair : eliminatedGenes) {
          pair.second /= fitnessSum;
        }
      }
    }

    // constexpr double NChooseK = 1;
    const double NChooseK = ::heu::NchooseK<double>(prevPopSize, K);
    // assert(!std::isnan(NChooseK));

    double newProbSum = 0;
    for (auto& pair : eliminatedGenes) {
      pair.second =
          NChooseK * std::pow(pair.second, K) * std::pow(1.0 - pair.second, prevPopSize - K);
      newProbSum += pair.second;
    }

    // assert(false);
    while (eliminatedGenes.size() > eliminateNum) {
      double r = randD(0.0, newProbSum);
      for (auto it = eliminatedGenes.begin(); it != eliminatedGenes.end(); ++it) {
        r -= it->second;
        if (r <= 0) {
          newProbSum -= it->second;
          eliminatedGenes.erase(it);
          break;
        }
      }

      if (r > 0) {
        newProbSum -= eliminatedGenes.back().second;
        eliminatedGenes.pop_back();
      }
    }

    // erase all eliminated candidates from the linked list
    for (auto& pair : eliminatedGenes) {
      static_cast<this_t*>(this)->_population.erase(pair.first);
    }

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(
        static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
  }
};

}  // namespace internal
}  // namespace heu

#endif  // HEU_SOGASELECTORS_HPP