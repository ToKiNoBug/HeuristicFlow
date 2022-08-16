/*
 Copyright © 2021-2022  TokiNoBug
This file is part of HeuristicFlow.

    HeuristicFlow is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HeuristicFlow is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HeuristicFlow.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef HEU_SOGASELECTORS_HPP
#define HEU_SOGASELECTORS_HPP

#include <unordered_map>

#include "GABase.hpp"

namespace heu {
namespace internal {

/**
 * \ingroup HEU_GENETIC
 * \class SOGASelector
 * \brief This class implement 10 different selection operators for SOGA. All kinds of selection
 * methods are listed in the enumeration type ``SelectMethod`.
 *
 * \tparam sm The selection method.
 */

template <SelectMethod sm>
class SOGASelector {
  // static_assert(false);
};

template <SelectMethod sm>
class SOGAInheriter;

template <>
class SOGASelector<SelectMethod::Truncation> {
  friend class SOGAInheriter<SelectMethod::Truncation>;

 protected:
  template <class this_t>
  void __impl___impl_select() noexcept {
    using Gene_t = typename this_t::Gene_t;

    using GeneIt_t = typename this_t::GeneIt_t;

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
  friend class SOGAInheriter<SelectMethod::RouletteWheel>;

 public:
  inline constexpr SelectMethod selectMethod() const noexcept {
    return SelectMethod::RouletteWheel;
  }

 protected:
  template <class this_t>
  void __impl___impl_select() noexcept {
    using GeneIt_t = typename this_t::GeneIt_t;

    /*
     * In this function, the linked list `eliminateSpace` stores all candidates to be eliminated.
     * At first it has all candidates, and as the loop running, we erase a candidate from
     * `eliminateSpace` where stores candidates waiting to be eliminated, symbolizing this
     * candidate is selected.
     *
     * When the size of `eliminateSpace` equals to eliminate num, all candidates inside
     * this hashmap will be erased from the linked list `_population`
     */

    std::list<std::pair<GeneIt_t, double>>
        eliminateSpace;  // Use the iterator as first and processed fitness as second

    // number of candidates that need to be eliminated.
    const int eliminateNum = int(static_cast<this_t*>(this)->_population.size() -
                                 static_cast<this_t*>(this)->_option.populationSize);

    // A cache to store the best fitness value in the previous generator cause the best gene may be
    // eliminated
    const double previousBestFitness = static_cast<this_t*>(this)->bestFitness();

    if (eliminateNum <= 0) {
      static_cast<this_t*>(this)->updateFailTimesAndBestGene(
          static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
      return;
    }
    // eliminateSpace.reserve(static_cast<this_t*>(this)->_population.size());
    double minFitness = pinfD;
    // fill eliminateSpace and compute the min fitness
    for (GeneIt_t it = static_cast<this_t*>(this)->_population.begin();
         it != static_cast<this_t*>(this)->_population.end(); ++it) {
      // If fitness option is FITNESS_LESS_BETTER, take the inverse value
      if constexpr (this_t::FitnessOpt == FitnessOption::FITNESS_GREATER_BETTER) {
        eliminateSpace.emplace_back(std::make_pair(it, it->fitness));
      } else {
        eliminateSpace.emplace_back(std::make_pair(it, -it->fitness));
      }
      minFitness = std::min(minFitness, eliminateSpace.back().second);
    }

    double processedFitnessSum = 0;
    for (auto& pair : eliminateSpace) {
      pair.second -= minFitness;
      processedFitnessSum += pair.second;
    }

    // erase all selected candidates from eliminateSpace

    while (int(eliminateSpace.size()) > eliminateNum) {
      // erase a solution from eliminateSpace by its fitness
      if (processedFitnessSum > 0) {
        double r = randD(0, processedFitnessSum);
        for (auto it = eliminateSpace.begin(); it != eliminateSpace.end(); ++it) {
          r -= it->second;
          if (r <= 0) {
            processedFitnessSum -= it->second;
            eliminateSpace.erase(it);
            break;
          }
        }

        if (r > 0) {  // unlikely
          processedFitnessSum -= eliminateSpace.back().second;
          eliminateSpace.pop_back();
        }

      } else {
        // select stochastically
        const int sizeBeforeErasement = int(eliminateSpace.size());
        for (auto it = eliminateSpace.begin(); it != eliminateSpace.end(); ++it) {
          if (randIdx(sizeBeforeErasement) == 0) {
            eliminateSpace.erase(it);
            break;
          }
        }
      }
    }

    static_cast<this_t*>(this)->applyErasement(eliminateSpace);

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(
        static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
  }
};

template <>
class SOGASelector<SelectMethod::Tournament> {
  friend class SOGAInheriter<SelectMethod::Tournament>;

 public:
  inline constexpr SelectMethod selectMethod() const noexcept { return SelectMethod::Tournament; }

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
    using Gene_t = typename this_t::Gene_t;
    using GeneIt_t = typename this_t::GeneIt_t;
    ;

    {
      const bool tournament_size_should_be_less_than_the_population_size =
          _tournamentSize < int(static_cast<this_t*>(this)->_option.populationSize);
      assert(tournament_size_should_be_less_than_the_population_size);
    }

    const int prevPopSize = int(static_cast<this_t*>(this)->_population.size());
    const double previousBestFitness = static_cast<this_t*>(this)->_bestGene->fitness;

    if (prevPopSize <= int(static_cast<this_t*>(this)->_option.populationSize)) {
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
    for (int playTimes = 0; playTimes < int(static_cast<this_t*>(this)->_option.populationSize);
         playTimes++) {
      //   find the best gene inside the tournament
      Gene_t* bestGenePtr = tournamentSpace.front();
      for (int idx = 1; idx < _tournamentSize; idx++) {
        if (this_t::isBetter(tournamentSpace[idx]->fitness, bestGenePtr->fitness)) {
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
      if (this_t::isBetter(it->fitness, curBestGene->fitness)) {
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
  friend class SOGAInheriter<SelectMethod::MonteCarlo>;

 public:
  inline constexpr SelectMethod selectMethod() const noexcept { return SelectMethod::MonteCarlo; }

 protected:
  template <class this_t>
  void __impl___impl_select() noexcept {
    using GeneIt_t = typename this_t::GeneIt_t;
    ;
    const double previousBestFitness = static_cast<this_t*>(this)->_bestGene->fitness;

    const int popSizeBeforeSelection = int(static_cast<this_t*>(this)->_population.size());

    if (popSizeBeforeSelection <= int(static_cast<this_t*>(this)->_option.populationSize)) {
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
        popSizeBeforeSelection - int(static_cast<this_t*>(this)->_option.populationSize);
    for (int idx = 0; idx < eliminateNum; idx++) {
      static_cast<this_t*>(this)->_population.erase(iterators[idx]);
    }

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(
        static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
  }
};

template <>
class SOGASelector<SelectMethod::Probability> {
  friend class SOGAInheriter<SelectMethod::Probability>;

 public:
  inline constexpr SelectMethod selectMethod() const noexcept { return SelectMethod::Probability; }

 protected:
  template <class this_t>
  void __impl___impl_select() noexcept {
    using GeneIt_t = typename this_t::GeneIt_t;

    std::list<std::pair<GeneIt_t, double>>
        eliminateSpace;  // Use the iterator as first and processed fitness as second
    const int prevPopSize = int(static_cast<this_t*>(this)->_population.size());
    const int K = int(static_cast<this_t*>(this)->_option.populationSize);
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
          eliminateSpace.emplace_back(std::make_pair(it, it->fitness));
        } else {
          eliminateSpace.emplace_back(std::make_pair(it, -it->fitness));
        }

        minFitness = std::min(eliminateSpace.back().second, minFitness);
      }

      //
      double fitnessSum = 0;
      for (auto& pair : eliminateSpace) {
        pair.second -= minFitness;
        fitnessSum += pair.second;
      }

      //
      if (fitnessSum <= 0) {
        for (auto& pair : eliminateSpace) {
          pair.second = 1.0 / prevPopSize;
        }
      } else {
        for (auto& pair : eliminateSpace) {
          pair.second /= fitnessSum;
        }
      }
    }

    constexpr double NChooseK = 1e100;
    // const double NChooseK = ::heu::NchooseK<double>(prevPopSize, K);
    //  assert(!std::isnan(NChooseK));

    double newProbSum = 0;
    for (auto& pair : eliminateSpace) {
      pair.second =
          NChooseK * std::pow(pair.second, K) * std::pow(1.0 - pair.second, prevPopSize - K);
      newProbSum += pair.second;
    }

    // assert(false);
    while (int(eliminateSpace.size()) > eliminateNum) {
      double r = randD(0.0, newProbSum);
      for (auto it = eliminateSpace.begin(); it != eliminateSpace.end(); ++it) {
        r -= it->second;
        if (r <= 0) {
          newProbSum -= it->second;
          eliminateSpace.erase(it);
          break;
        }
      }

      if (r > 0) {
        newProbSum -= eliminateSpace.back().second;
        eliminateSpace.pop_back();
      }
    }

    static_cast<this_t*>(this)->applyErasement(eliminateSpace);

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(
        static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
  }
};

template <>
class SOGASelector<SelectMethod::LinearRank> {
  friend class SOGAInheriter<SelectMethod::LinearRank>;

 public:
  SOGASelector() {
    _linearSelectBestProbability = 0.9;
    _linearSelectWrostProbability = 0.1;
  }

  inline constexpr SelectMethod selectMethod() const noexcept { return SelectMethod::LinearRank; }

  inline double linearSelectBestProbability() const noexcept {
    return _linearSelectBestProbability;
  }

  inline double linearSelectWrostProbability() const noexcept {
    return _linearSelectWrostProbability;
  }

  inline void setLinearSelectProbability(double worstProbabilty, double bestProbability) noexcept {
    const bool worst_gene_select_probability_should_be_less_than_the_best_gene_select_probability =
        worstProbabilty < bestProbability;
    assert(worst_gene_select_probability_should_be_less_than_the_best_gene_select_probability);

    const bool probability_should_be_greater_or_equal_to_0 =
        worstProbabilty >= 0 && bestProbability >= 0;

    assert(probability_should_be_greater_or_equal_to_0);

    const bool probability_should_be_less_or_equal_to_1 =
        worstProbabilty <= 1 && bestProbability <= 1;
    assert(probability_should_be_less_or_equal_to_1);

    _linearSelectBestProbability = bestProbability;
    _linearSelectWrostProbability = worstProbabilty;
  }

 protected:
  double _linearSelectWrostProbability;
  double _linearSelectBestProbability;

  template <class this_t>
  void __impl___impl_select() noexcept {
    using GeneIt_t = typename this_t::GeneIt_t;

    const int popSizeBeforeSelect = int(static_cast<this_t*>(this)->_population.size());
    const int K = int(static_cast<this_t*>(this)->_option.populationSize);
    const double previousBestFitness = static_cast<this_t*>(this)->bestFitness();

    if (popSizeBeforeSelect <= K) {
      static_cast<this_t*>(this)->updateFailTimesAndBestGene(
          static_cast<this_t*>(this)->findCurrentBestGene());
      return;
    }

    std::vector<GeneIt_t> sortSpace(0);
    sortSpace.reserve(popSizeBeforeSelect);
    for (GeneIt_t it = static_cast<this_t*>(this)->_population.begin();
         it != static_cast<this_t*>(this)->_population.end(); ++it) {
      sortSpace.emplace_back(it);
    }

    std::sort(sortSpace.begin(), sortSpace.end(), this_t::GeneItCompareFun);

    const double nNegitive = _linearSelectWrostProbability * popSizeBeforeSelect;
    const double nPositive = _linearSelectBestProbability * popSizeBeforeSelect;

    std::list<std::pair<GeneIt_t, double>> eliminateSpace;
    for (int idx = 0; idx < int(sortSpace.size()); idx++) {
      eliminateSpace.emplace_back(std::make_pair(
          sortSpace[idx],

          (nNegitive +
           (nPositive - nNegitive) * (popSizeBeforeSelect - idx - 1) / (popSizeBeforeSelect - 1)) /
              popSizeBeforeSelect

          ));
    }

    double rMax = 1.0;
    while (int(eliminateSpace.size()) > popSizeBeforeSelect - K) {
      double r = randD(0, rMax);
      for (auto it = eliminateSpace.begin(); it != eliminateSpace.end(); ++it) {
        r -= it->second;
        if (r <= 0) {
          rMax -= it->second;
          eliminateSpace.erase(it);
          break;
        }
      }

      if (r > 0) {
        rMax -= eliminateSpace.back().second;
        eliminateSpace.pop_back();
      }
    }

    static_cast<this_t*>(this)->applyErasement(eliminateSpace);

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(
        static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
  }
};

template <>
class SOGASelector<SelectMethod::ExponentialRank> {
  friend class SOGAInheriter<SelectMethod::ExponentialRank>;

 public:
  SOGASelector() { _exponetialSelectBase = 0.8; }

  inline constexpr SelectMethod selectMethod() const noexcept {
    return SelectMethod::ExponentialRank;
  }

  inline double exponetialSelectBase() const noexcept { return _exponetialSelectBase; }

  inline void setExponetialSelectBase(double _c) noexcept {
    const bool the_base_number_for_exponential_rank_selection_should_be_greater_or_equal_to_0 =
        _c >= 0;
    assert(the_base_number_for_exponential_rank_selection_should_be_greater_or_equal_to_0);

    const bool the_base_number_for_exponential_rank_selection_should_be_less_than_1 = _c < 1;
    assert(the_base_number_for_exponential_rank_selection_should_be_less_than_1);

    _exponetialSelectBase = _c;
  }

 protected:
  double _exponetialSelectBase;

  template <class this_t>
  void __impl___impl_select() noexcept {
    using GeneIt_t = typename this_t::GeneIt_t;
    using GeneIt_t = typename this_t::GeneIt_t;

    const int popSizeBeforeSelect = int(static_cast<this_t*>(this)->_population.size());
    const int K = int(static_cast<this_t*>(this)->_option.populationSize);
    const double previousBestFitness = static_cast<this_t*>(this)->bestFitness();

    if (popSizeBeforeSelect <= K) {
      static_cast<this_t*>(this)->updateFailTimesAndBestGene(
          static_cast<this_t*>(this)->findCurrentBestGene());
      return;
    }

    std::vector<GeneIt_t> sortSpace(0);
    sortSpace.reserve(popSizeBeforeSelect);
    for (GeneIt_t it = static_cast<this_t*>(this)->_population.begin();
         it != static_cast<this_t*>(this)->_population.end(); ++it) {
      sortSpace.emplace_back(it);
    }

    std::sort(sortSpace.begin(), sortSpace.end(), this_t::GeneItCompareFun);

    const double c_minus_1_div_c_pow_N_minus_1 =
        (_exponetialSelectBase - 1) / (std::pow(_exponetialSelectBase, popSizeBeforeSelect) - 1);

    std::list<std::pair<GeneIt_t, double>> eliminateSpace;
    for (int idx = 0; idx < int(sortSpace.size()); idx++) {
      eliminateSpace.emplace_back(std::make_pair(
          sortSpace[idx],

          c_minus_1_div_c_pow_N_minus_1 * std::pow(_exponetialSelectBase, idx)

              ));
    }

    double rMax = 1.0;
    while (int(eliminateSpace.size()) > popSizeBeforeSelect - K) {
      double r = randD(0, rMax);
      for (auto it = eliminateSpace.begin(); it != eliminateSpace.end(); ++it) {
        r -= it->second;
        if (r <= 0) {
          rMax -= it->second;
          eliminateSpace.erase(it);
          break;
        }
      }

      if (r > 0) {
        rMax -= eliminateSpace.back().second;
        eliminateSpace.pop_back();
      }
    }

    static_cast<this_t*>(this)->applyErasement(eliminateSpace);

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(
        static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
  }
};

template <>
class SOGASelector<SelectMethod::Boltzmann> {
  friend class SOGAInheriter<SelectMethod::Boltzmann>;

 public:
  inline constexpr SelectMethod selectMethod() const noexcept { return SelectMethod::Boltzmann; }

  inline double boltzmannSelectStrength() noexcept { return _boltzmannSelectStrength; }

  inline void setBoltzmannSelectStrength(double b) noexcept { _boltzmannSelectStrength = b; }

 protected:
  double _boltzmannSelectStrength;

  template <class this_t>
  void __impl___impl_select() noexcept {
    using GeneIt_t = typename this_t::GeneIt_t;

    // static_cast<this_t*>(this)

    const double previousBestFitness = static_cast<this_t*>(this)->_bestGene->fitness;

    const int popSizeBeforeSelect = int(static_cast<this_t*>(this)->_population.size());

    const int eliminateNum =
        popSizeBeforeSelect - int(static_cast<this_t*>(this)->_option.populationSize);

    if (eliminateNum <= 0) {
      static_cast<this_t*>(this)->updateFailTimesAndBestGene(
          static_cast<this_t*>(this)->findCurrentBestGene());
      return;
    }

    std::list<std::pair<GeneIt_t, double>> eliminateSpace;
    // double Z = 0;    //Z is erased because probability normalization is not guaranteed.
    for (GeneIt_t it = static_cast<this_t*>(this)->_population.begin();
         it != static_cast<this_t*>(this)->_population.end(); ++it) {
      eliminateSpace.emplace_back(it, it->fitness);
      // Z += std::exp(it->fitness);
    }

    double rMax = 0;  // the sum of "probability"
    for (auto& pair : eliminateSpace) {
      pair.second = std::exp(_boltzmannSelectStrength * pair.second);
      rMax += pair.second;  //  the sum of probability may not be 1.
    }

    while (int(eliminateSpace.size()) > eliminateNum) {
      double r = randD(0, rMax);

      for (auto it = eliminateSpace.begin(); it != eliminateSpace.end(); ++it) {
        r -= it->second;
        if (r <= 0) {
          rMax -= it->second;
          eliminateSpace.erase(it);
          break;
        }
      }

      if (r > 0) {
        rMax -= eliminateSpace.back().second;
        eliminateSpace.pop_back();
      }
    }

    static_cast<this_t*>(this)->applyErasement(eliminateSpace);

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(
        static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
  }
};

template <>
class SOGASelector<SelectMethod::StochasticUniversal> {
  friend class SOGAInheriter<SelectMethod::StochasticUniversal>;

 protected:
  template <class this_t>
  void __impl___impl_select() noexcept {
    using GeneIt_t = typename this_t::GeneIt_t;

    const double previousBestFitness = static_cast<this_t*>(this)->_bestGene->fitness;

    const int popSizeBeforeSelect = int(static_cast<this_t*>(this)->_population.size());

    const int eliminateNum =
        popSizeBeforeSelect - int(static_cast<this_t*>(this)->_option.populationSize);

    if (eliminateNum <= 0) {
      static_cast<this_t*>(this)->updateFailTimesAndBestGene(
          static_cast<this_t*>(this)->findCurrentBestGene());
      return;
    }

    std::list<std::pair<GeneIt_t, double>> eliminateSpace;

    double minFitness = pinfD;

    for (GeneIt_t it = static_cast<this_t*>(this)->_population.begin();
         it != static_cast<this_t*>(this)->_population.end(); ++it) {
      // If fitness option is FITNESS_LESS_BETTER, take the inverse value
      if constexpr (this_t::FitnessOpt == FitnessOption::FITNESS_GREATER_BETTER) {
        eliminateSpace.emplace_back(std::make_pair(it, it->fitness));
      } else {
        eliminateSpace.emplace_back(std::make_pair(it, -it->fitness));
      }

      minFitness = std::min(minFitness, eliminateSpace.back().second);
    }

    double fitnessSum = 0;
    for (auto& pair : eliminateSpace) {
      pair.second -= minFitness;

      fitnessSum += pair.second;
    }

    if (fitnessSum <= 0) {  //  unlikely
      for (auto it = eliminateSpace.begin(); it != eliminateSpace.end();) {
        if (randD() * static_cast<this_t*>(this)->_option.populationSize <= popSizeBeforeSelect) {
          it = eliminateSpace.erase(it);
        } else {
          ++it;
        }
      }
    } else {
      const double pointerInterval =
          fitnessSum / static_cast<this_t*>(this)->_option.populationSize;
      /*
      const double begPosition = randD(0, pointerInterval);

      double passedFitnessSum = begPositon;
      auto it = eliminateSpace.begin();
      for (int pointerIdx = 0; pointerIdx < static_cast<this_t*>(this)->_option.populationSize;
           pointerIdx++) {
        if (it == eliminateSpace.end()) break;  //  Theoritically this won't happen

        while (passedFitnessSum < (pointerIdx + 1) * pointerInterval) {
          passedFitnessSum += it->second;
          ++it;
        }
        it = eliminateSpace.erase(it);
      */

      // While paving throught eliminateSpace, accumulate the fitness. Everytime when the sum is
      // greater of equal to the interval, erase the current element (selected)
      double fitnessCounter = randD(0, pointerInterval);
      for (auto it = eliminateSpace.begin(); it != eliminateSpace.end();) {
        fitnessCounter += it->second;

        if (fitnessCounter >= pointerInterval) {
          fitnessCounter -= pointerInterval;
          it = eliminateSpace.erase(it);
        } else {
          ++it;
        }
      }

      while (int(eliminateSpace.size()) > eliminateNum) {
        //  This rarely happens, but we must ensure
        //  the size of population keep unchanged in each generation
        eliminateSpace.pop_back();
      }
    }

    static_cast<this_t*>(this)->applyErasement(eliminateSpace);

    static_cast<this_t*>(this)->updateFailTimesAndBestGene(
        static_cast<this_t*>(this)->findCurrentBestGene(), previousBestFitness);
  }
};

template <>
class SOGASelector<SelectMethod::EliteReserved> {
  friend class SOGAInheriter<SelectMethod::EliteReserved>;

 public:
  SOGASelector() { _eliteNum = 1; }

  inline constexpr SelectMethod selectMethod() const noexcept {
    return SelectMethod::EliteReserved;
  }

  inline int eliteNum() const noexcept { return _eliteNum; }

  inline void setEliteNum(int eliteNum) noexcept {
    assert(eliteNum >= 1);
    _eliteNum = eliteNum;
  }

 protected:
  int _eliteNum;

  template <class this_t>
  void __impl___impl_select() noexcept {
    using GeneIt_t = typename this_t::GeneIt_t;
    const int popSizeBeforeSelect = int(static_cast<this_t*>(this)->_population.size());
    const int popSizeAfterSelect = int(static_cast<this_t*>(this)->_option.populationSize);
    const int eliminateNum = popSizeBeforeSelect - popSizeAfterSelect;

    if (eliminateNum <= 0) {
      static_cast<this_t*>(this)->updateFailTimesAndBestGene(
          static_cast<this_t*>(this)->findCurrentBestGene());
      return;
    }

    std::list<std::pair<GeneIt_t, double>> eliminateSpace;
    double minFitness = pinfD;
    {
      std::vector<GeneIt_t> sortSpace(0);

      sortSpace.reserve(popSizeBeforeSelect);

      for (GeneIt_t it = static_cast<this_t*>(this)->_population.begin();
           it != static_cast<this_t*>(this)->_population.end(); ++it) {
        sortSpace.emplace_back(it);
      }

      std::sort(sortSpace.begin(), sortSpace.end(), this_t::GeneItCompareFun);

      for (int idx = _eliteNum; idx < int(sortSpace.size()); idx++) {
        if constexpr (this_t::FitnessOpt == FitnessOption::FITNESS_GREATER_BETTER)
          eliminateSpace.emplace_back(sortSpace[idx], sortSpace[idx]->fitness);
        else
          eliminateSpace.emplace_back(sortSpace[idx], -sortSpace[idx]->fitness);

        minFitness = std::min(minFitness, eliminateSpace.back().second);
      }
    }

    double rMax = 0;

    for (auto& pair : eliminateSpace) {
      pair.second -= minFitness;
      rMax += pair.second;
    }

    if (rMax <= 0) {  //  erase randomly
      while (int(eliminateSpace.size()) > eliminateNum) {
        const int curEliminateSpaceSize = int(eliminateSpace.size());
        for (auto it = eliminateSpace.begin(); it != eliminateSpace.end(); ++it) {
          if (randD() * curEliminateSpaceSize <= 1) {
            eliminateSpace.erase(it);
            break;
          }
        }
      }
    } else {
      while (int(eliminateSpace.size()) > eliminateNum) {
        double r = randD(0, rMax);

        for (auto it = eliminateSpace.begin(); it != eliminateSpace.end(); ++it) {
          r -= it->second;
          if (r <= 0) {
            rMax -= it->second;
            eliminateSpace.erase(it);
            break;
          }
        }
        if (r > 0) {  //  unlikely
          rMax -= eliminateSpace.back().second;
          eliminateSpace.pop_back();
        }
      }
    }

    static_cast<this_t*>(this)->applyErasement(eliminateSpace);
    static_cast<this_t*>(this)->updateFailTimesAndBestGene(
        static_cast<this_t*>(this)->findCurrentBestGene());
  }
};

template <SelectMethod sm>
class SOGAInheriter : public SOGAInheriter<SelectMethod((int)(sm) + 1)> {
  using next = SOGAInheriter<SelectMethod((int)(sm) + 1)>;

 public:
  struct type : public SOGASelector<sm>, public next::type {};

  template <class this_t>
  inline static void __impl___impl___impl_select(type* solver, const SelectMethod _sm) noexcept {
    if constexpr (sm + 1 >= SelectMethod::RunTimeSelectMethod) {
      static_cast<SOGASelector<sm>*>(solver)->template __impl___impl_select<this_t>();
    } else {
      if (_sm == sm) {
        static_cast<SOGASelector<sm>*>(solver)->template __impl___impl_select<this_t>();
      } else {
        next::template __impl___impl___impl_select<this_t>(solver, _sm);
      }
    }
  }
};

template <>
class SOGAInheriter<SelectMethod::RunTimeSelectMethod> {
 public:
  struct type {};
};

template <>
class SOGASelector<SelectMethod::RunTimeSelectMethod>
    : public SOGAInheriter<SelectMethod(0)>::type {
 public:
  SOGASelector() { _selectMethod = SelectMethod::EliteReserved; }
  using someType = SOGAInheriter<SelectMethod(0)>;
  inline SelectMethod selectMethod() const noexcept { return _selectMethod; }

  inline void setSelectMethod(SelectMethod _sm) noexcept {
    switch (_sm) {
      case SelectMethod::Boltzmann:
      case SelectMethod::EliteReserved:
      case SelectMethod::ExponentialRank:
      case SelectMethod::LinearRank:
      case SelectMethod::MonteCarlo:
      case SelectMethod::Probability:
      case SelectMethod::RouletteWheel:
      case SelectMethod::StochasticUniversal:
      case SelectMethod::Tournament:
      case SelectMethod::Truncation:
        _selectMethod = _sm;
        break;
      default:
        const bool invalid_select_method = false;
        assert(invalid_select_method);
        break;
    }
  }

 protected:
  SelectMethod _selectMethod;

  template <class this_t>
  inline void __impl___impl_select() noexcept {
    SOGAInheriter<SelectMethod(0)>::template __impl___impl___impl_select<this_t>(this,
                                                                                 _selectMethod);
  }
};

}  // namespace internal
}  // namespace heu

#endif  // HEU_SOGASELECTORS_HPP
