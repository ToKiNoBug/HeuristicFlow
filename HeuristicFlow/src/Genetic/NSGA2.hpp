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

#ifndef HEU_NSGA2BASE_HPP
#define HEU_NSGA2BASE_HPP

#include "InternalHeaderCheck.h"
#include "NSGABase.hpp"

namespace heu {

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class NSGA2
 * \brief NSGA2 MOGA solver. Suitable for not too many objectives.
 *
 * This class implemented the NSGA-II algorithm.\n
 * NSGA-II is a template for multi-objective genetic algorithm based on Pareto optimality. It selects by nondomainance
 * sorting.
 *
 * \sa NSGA2::select for this special procedure.
 *
 * @tparam Var_t Type of decisition variable.
 * @tparam ObjNum Numbers of objectives.
 * @tparam fOpt Whether greater fitness value means better.
 * @tparam rOPt Whether the solver records fitness changelog.
 * @tparam Args_t Type of other parameters.
 * @tparam _iFun_ Compile-time iFun, use nullptr for runtime
 * @tparam _fFun_ Compile-time fFun, use nullptr for runtime
 * @tparam _cFun_ Compile-time cFun, use nullptr for runtime
 * @tparam _mFun_ Compile-time mFun, use nullptr for runtime
 *
 * \sa GAOption for ga running parameters
 * \sa SOGA for APIs that all genetic solvers have.
 * \sa NSGA3 for many objective problems
 *
 * ## APIs that MOGA solvers have:
 * - `void paretoFront(std::vector<Fitness_t>& front) const` get pareto front of fitness values.
 * - `void paretoFront(std::vector<std::pair<const Var_t*, const Fitness_t*>>& front) const` get pareto front of
 * decision variables and fitness values.
 * - `const std::unordered_set<const Gene*>& pfGenes() const` returns a const reference to the PF hash set.
 * - `size_t objectiveNum() const` get the number of objectives.
 *
 * ## APIs that MOGA with dynamic objective numbers have:
 * - `void setObjectiveNum(int _objNum)` set the objective number.
 *
 */
template <typename Var_t, int ObjNum, FitnessOption fOpt = FITNESS_LESS_BETTER, RecordOption rOpt = DONT_RECORD_FITNESS,
          class Args_t = void,
          typename internal::GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::initializeFun _iFun_ = nullptr,
          typename internal::GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::fitnessFun _fFun_ = nullptr,
          typename internal::GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::crossoverFun _cFun_ = nullptr,
          typename internal::GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::mutateFun _mFun_ = nullptr>
class NSGA2 : public internal::NSGABase<Var_t, ObjNum, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
  using Base_t = internal::NSGABase<Var_t, ObjNum, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

 private:
 public:
  NSGA2(){};
  ~NSGA2(){};
  HEU_MAKE_NSGABASE_TYPES(Base_t)
  friend class internal::GABase<Var_t, Fitness_t, DONT_RECORD_FITNESS, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

  /**
   * \brief This struct is used to store nondomainance sorting-related informations in select operation. A pointer to
   * this struct, infoUnit2*, has access to every information to a gene.
   *
   * There's several sorting by according to different attribute, thus a vector of infoUnit2* is widely used when
   * sorting.
   *
   */
  struct infoUnit2 : public infoUnitBase_t {
   public:
    /**
     * \brief Congestion value of this gene.
     *
     */
    double congestion;
  };

  /**
   * \brief Run the solver.
   *
   */
  inline void run() { this->template __impl_run<NSGA2>(); }

 protected:
  /**
   * \brief Sort according to infoUnitBase::domainedByNum in an acsending order. This sorting step will seperate the
   * whole population into several pareto layers.
   *
   * \param A
   * \param B
   * \return true A's congestion is greater than B
   * \return false A's congestion is not greater than B
   */
  static bool compareByCongestion(const infoUnitBase_t *A, const infoUnitBase_t *B) {
    if (A == B) return false;
    return (static_cast<const infoUnit2 *>(A)->congestion) > (static_cast<const infoUnit2 *>(B)->congestion);
  }

  /**
   * \brief Template parameter objIdx is used to indicate how to compare 2 const infoUnit *.
   *
   * \tparam objIdx Indicate which objective is compared when sorting.
   * \param A
   * \param B
   * \return true objIdx's fitness of A is greater than B.
   * \return false objIdx's fitness of A is not greater than B.
   *
   * \note This procedure actually has no relationship with template parameter fOpt, because sorting as fitness is
   * simply used to compute congestion.
   */
  template <int64_t objIdx>
  static bool compareByFitness(const infoUnitBase_t *A, const infoUnitBase_t *B) {
    static_assert(objIdx >= 0, "Invalid comparison flag");
    if (A == B) return false;
    /// compare by fitness on single objective
    return A->fitnessCache[objIdx] < B->fitnessCache[objIdx];
  }

  /**
   * \brief Non-dominated sorting
   *
   * This is the core of NSGA2.
   *
   * Simply put, the whole population is divided into several non-dominated layers, relatively front layers will be
   * selected. Usually there will be a layer that can't be selected entirely, we must select some and eliminate the
   * rest. To achieve this, we compute congestion so that genes will greater congestion are more likely to be selected.
   * That's how NSGA2 maintain variety in the PF.
   *
   * Besides, as the algorithm running, it's common to see that the whole population after selection are all members in
   * PF. Don't worry, NSGA2 works well with that condition.
   *
   * In detail, this function applies non-dominated sorting in following steps:
   * 1. Make a `std::vector` of `infoUnit2` named `pop` that each element corresponds to an element in
   * population(`std::list<Gene>`).
   * 2. Fill `this->sortSpace` (std::list<infoUnitBase*) with pointer to each member in `pop`.
   * 3. Compute `infoUnit::domainedByNum` for `pop`.
   * 4. Sort elements in `this->sortSpace` according to `infoUnitBase::domainedByNum`.
   * 5. Divide the whole population into several nondominated layers. This function is implemented in `NSGABase`, and
   * the result is stored in `this->pfLayers`. For better performance, each `std::vector<infoUnit *>` will reserve space
   * for the whole population. Here we use vector since sorting might happen in a single layer.
   * 6. Make a `std::unordered_set` of `infoUnit2 *` named `selected` to store all selected genes. Insert the Pareto
   * frontier and each layers into `selected` until a layer can't be completedly inserted(I name this layer *K*).
   * Congestion are used to select several genes in *K*. Besides, there is a possible that a layer is completed inserted
   * and size of `selected` equals to user assigned population size. In that condition, jump directly to step 10.
   * 7. Sort the whole population according to every objective values to compute the congestion.
   * 8. Sort elements in *K* by its congestion in descend order. Elements with greater congestion have less index.
   * 9. Inserts each element in *K* into `selected` (by index order) until size of `selected` equals to user assigned
   * population size(`GABase::_option.populationSize`).
   * 10. Go through the whole population and erase if a gene doesn't exist in `selected`.


   * Every sorting steps use `std::sort` and different comparision function are used to sort according to different
   * attributes of a gene.
   *
   */
  void __impl_select() {
    using cmpFun_t = bool (*)(const infoUnitBase_t *, const infoUnitBase_t *);
    static const size_t objCapacity = (ObjNum == Eigen::Dynamic) ? (HEU_MAX_RUNTIME_OBJNUM) : ObjNum;
    static const std::array<cmpFun_t, objCapacity> fitnessCmpFuns = expand<0, objCapacity - 1>();

    const size_t popSizeBefore = this->_population.size();
    std::vector<infoUnit2> pop;
    pop.clear();
    pop.reserve(popSizeBefore);
    this->sortSpace.resize(popSizeBefore);

    for (auto it = this->_population.begin(); it != this->_population.end(); ++it) {
      pop.emplace_back();
      pop.back().iterator = it;
      pop.back().fitnessCache = it->_Fitness;
      pop.back().congestion = 0;
    }

    // make sortspace
    for (size_t i = 0; i < popSizeBefore; i++) {
      this->sortSpace[i] = pop.data() + i;
    }

    this->calculateDominatedNum();

    this->divideLayers();

    const size_t PFSize = this->pfLayers.front().size();
    if (PFSize <= this->_option.populationSize)
      this->updatePF((const infoUnitBase_t **)(this->pfLayers.front().data()), this->pfLayers.front().size());

    std::unordered_set<infoUnit2 *> selected;
    selected.reserve(this->_option.populationSize);
    bool needCongestion = true;
    while (true) {
      // don't need to calculate congestion
      if (selected.size() == this->_option.populationSize) {
        needCongestion = false;
        break;
      }
      // need to calculate congestion
      if (selected.size() + this->pfLayers.front().size() > this->_option.populationSize) {
        needCongestion = true;
        break;
      }
      // emplace every element of this layer into selected
      for (const auto i : this->pfLayers.front()) {
        selected.emplace(static_cast<infoUnit2 *>(i));
      }
      this->pfLayers.pop_front();
    }
    // calculate congestion
    if (needCongestion) {
      for (size_t objIdx = 0; objIdx < this->objectiveNum(); objIdx++) {
        std::sort((infoUnit2 **)(this->sortSpace.data()),
                  (infoUnit2 **)(this->sortSpace.data() + this->sortSpace.size()), fitnessCmpFuns[objIdx]);

        const double scale =
            std::abs(this->sortSpace.front()->fitnessCache[objIdx] - this->sortSpace.back()->fitnessCache[objIdx]) +
            1e-10;

        static_cast<infoUnit2 *>(this->sortSpace.front())->congestion = internal::pinfD;
        static_cast<infoUnit2 *>(this->sortSpace.back())->congestion = internal::pinfD;

        // calculate congestion on single object
        for (size_t idx = 1; idx < popSizeBefore - 1; idx++) {
          static_cast<infoUnit2 *>(this->sortSpace[idx])->congestion +=
              std::abs(this->sortSpace[idx - 1]->fitnessCache[objIdx] -
                       this->sortSpace[idx + 1]->fitnessCache[objIdx]) /
              scale;
        }
      }  // end sort on objIdx

      // sort by congestion in the undetermined set
      std::sort((infoUnit2 **)(this->pfLayers.front().data()),
                (infoUnit2 **)(this->pfLayers.front().data() + this->pfLayers.front().size()), compareByCongestion);

      size_t idx = 0;
      while (selected.size() < this->_option.populationSize) {
        selected.emplace(static_cast<infoUnit2 *>(this->pfLayers.front()[idx]));
        idx++;
      }

    }  // end applying congestion

    // erase unselected
    for (auto &i : pop) {
      if (selected.find(&i) == selected.end()) {
        this->_population.erase(i.iterator);
      }
    }

    if (PFSize > this->_option.populationSize) {
      std::vector<const infoUnitBase_t *> PF;
      PF.reserve(selected.size());
      for (auto i : selected) {
        PF.emplace_back(i);
      }
      this->updatePF(PF.data(), PF.size());
    }

    this->pfLayers.clear();
    this->sortSpace.clear();
  }

 private:
  // some template metaprogramming to make a function pointer array as below:
  // universialCompareFun<0>,universialCompareFun<1>,...,universialCompareFun<ObjNum-1>
  using fun_t = bool (*)(const infoUnitBase_t *, const infoUnitBase_t *);

  template <int64_t cur, int64_t max>
  struct expandStruct {
    static void expand(fun_t *dst) {
      *dst = compareByFitness<cur>;
      expandStruct<cur + 1, max>::expand(dst + 1);
    }
  };

  template <int64_t max>
  struct expandStruct<max, max> {
    static void expand(fun_t *dst) { *dst = compareByFitness<max>; }
  };

  template <int64_t beg, int64_t end>
  std::array<fun_t, end - beg + 1> expand() {
    std::array<fun_t, end - beg + 1> funs;
    expandStruct<beg, end>::expand(funs.data());
    return funs;
  }
};

}  // namespace heu

#endif  // HEU_NSGA2BASE_HPP
