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
#ifndef HEU_NSGA3ABSTRACT_HPP
#define HEU_NSGA3ABSTRACT_HPP

#include <unordered_map>
#include <unordered_set>

#include "InternalHeaderCheck.h"
#include "NSGABase.hpp"

namespace heu {

namespace internal {

/**
 * \ingroup HEU_GENETIC
 * \brief Internal base class for NSGA3.
 *
 * This class implements most part of NSGA3' selection precedure. Also it maintains reference points in a matrix.
 *
 * \tparam Var_t
 * \tparam ObjNum
 * \tparam rOpt
 * \tparam Args_t
 * \tparam _iFun_
 * \tparam _fFun_
 * \tparam _cFun_
 * \tparam _mFun_
 */
template <typename Var_t, int ObjNum, RecordOption rOpt, class Args_t,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::initializeFun _iFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::fitnessFun _fFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::crossoverFun _cFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::mutateFun _mFun_>
class NSGA3Abstract
    : public NSGABase<Var_t, ObjNum, FITNESS_LESS_BETTER, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
  using Base_t = NSGABase<Var_t, ObjNum, FITNESS_LESS_BETTER, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

 public:
  ~NSGA3Abstract() {}
  HEU_MAKE_NSGABASE_TYPES(Base_t)
  friend class internal::GABase<Var_t, Fitness_t, DONT_RECORD_FITNESS, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

  using RefPointIdx_t = size_t;

  /// Type of reference point(called RP in brief) matrix. Each coloumn is the coordinate of a RP.
  using RefMat_t = Eigen::Array<double, ObjNum, Eigen::Dynamic>;

  /**
   * \brief Get RPs
   *
   * \return const RefMat_t& Const reference to RPs
   */
  inline const RefMat_t& referencePoints() const { return referencePoses; }

  /**
   * \brief The struct to store all information about a gene used in NSGA3's selection.
   *
   */
  struct infoUnit3 : public infoUnitBase_t {
    /// The translated fitness. Normalized fitness is also stored in this.
    Fitness_t translatedFitness;
    /// The index of its closet RP
    size_t closestRefPoint;
    /// Distance to the closet RP
    double distance;
  };

 protected:
  /// RP matrix. Each coloumn is the coordinate of a RP.
  RefMat_t referencePoses;

  /**
   * \brief Generate a series of reference points using Das and Dennis’s method.
   *
   * \param dimN Number of objectives
   * \param precision Count of points, usually 3 to 5. Greater precision means much more RPs
   * \param dst Where to put the PF points temporarly.
   */
  void computeReferencePointPoses(const size_t dimN, const size_t precision, std::vector<Fitness_t>* dst) const {
    dst->clear();
    dst->reserve(NchooseK(dimN + precision - 1, precision));

    pri_startRP(dimN, precision, dst);
  }

  /**
   * \brief The core procedure of NSGA3.
   *
   * Simply put, NSGA3 uses reference points to maintain a diverse PF. It's simliar with NSGA2 but how to select
   * elements from the undetermined layer. The undetermined layer is the first non-dominated layer that can't be
   * selected entirely. In document of NSGA2, the undetermined layer is named *K*.
   *
   *
   * Detailed steps: (Steps same in NSGA2 won't be introduced again)
   * 1. Make a vector of `infoUnit3` and fill `this->sortSpace` (same as NSGA2)
   * 2. Compute dominated number of the population. (same as NSGA2)
   * 3. Divide the population into several layers. (same as NSGA2)
   * 4. Insert each layer into `selected` (`std::unordered_set<infoUnit3*>`) until a layer can't be inserted entirely.
   * (same as NSGA2).This undetermined layer is named `Fl`.
   * 5. Normalize procedure
   * 6. Associate procedure
   * 7. Niche preservation procedure.
   * 8. Erase all unselected genes.
   *
   * \sa NSGA2
   */
  void __impl_select() {
    // population size before selection.
    const size_t popSizeBef = this->_population.size();
    std::vector<infoUnit3> pop;
    pop.reserve(popSizeBef);
    for (auto it = this->_population.begin(); it != this->_population.end(); ++it) {
      pop.emplace_back();
      pop.back().iterator = it;
      pop.back().fitnessCache = it->_Fitness;
      pop.back().closestRefPoint = -1;
    }

    this->sortSpace.resize(popSizeBef);
    for (size_t i = 0; i < popSizeBef; i++) {
      this->sortSpace[i] = pop.data() + i;
    }

    this->calculateDominatedNum();
    this->divideLayers();

    const size_t PFSize = this->pfLayers.front().size();

    if (PFSize <= this->_option.populationSize)
      this->updatePF((const infoUnitBase_t**)this->pfLayers.front().data(), this->pfLayers.front().size());

    // hash set to store genes that will be selected
    std::unordered_set<infoUnit3*> selected;
    selected.reserve(this->_option.populationSize);

    // pointer to undertermined layer Fl
    std::vector<infoUnitBase_t*>* FlPtr = nullptr;

    // if need to use RP in this selection
    bool needRefPoint;

    while (true) {
      // Luckily we don't need to select part of some genes in Fl.
      if (selected.size() == this->_option.populationSize) {
        needRefPoint = false;
        break;
      }

      // We need to select part of members in Fl and eliminate the rest
      if (selected.size() + this->pfLayers.front().size() > this->_option.populationSize) {
        needRefPoint = true;
        FlPtr = &this->pfLayers.front();
        break;
      }

      // emplace the whole layer into selected
      for (infoUnitBase_t* i : this->pfLayers.front()) {
        selected.emplace(static_cast<infoUnit3*>(i));
      }
      this->pfLayers.pop_front();
    }

    if (needRefPoint) {
      // Normalize procedure
      std::unordered_multimap<RefPointIdx_t, infoUnit3*>
          Fl;  // Fl is stored in a multihash sothat we can find a gene by its related RP
      Fl.reserve(FlPtr->size());
      normalize(selected, *FlPtr);
      // use reference points' index(col idx) as key and genes related to this point as value.
      std::unordered_map<RefPointIdx_t, size_t> refPoints;
      refPoints.reserve(referencePoses.cols());
      for (int i = 0; i < referencePoses.cols(); i++) {
        refPoints.emplace(i, 0);
      }

      // Associate procedure
      associate(selected);
      associate(*FlPtr, &Fl);

      // niche preservation procedure.
      nichePreservation(&selected, &Fl, &refPoints);
    }

    // erase all unselected genes
    for (auto i : this->sortSpace) {
      if (selected.find((infoUnit3*)i) == selected.end()) {
        this->_population.erase(i->iterator);
      }
    }

    // update PF after selection if the next population is filled with PF(We don't hope the size of PF exceeds size of
    // population)
    if (PFSize > this->_option.populationSize) {
      std::vector<const infoUnitBase_t*> PF;
      PF.reserve(selected.size());
      for (auto i : selected) {
        PF.emplace_back(i);
      }
      this->updatePF(PF.data(), PF.size());
    }

  }  //  end selection

  /**
   * \brief Normalize procedure automatically normalize all objectives, since different objectives may have different
   * orders of magnitude.
   *
   * \param selected Genes that are already selected.
   * \param Fl Genes to be partly selected
   *
   * This procedure has following steps:
   * 1. Find a ideal point of (selected∪Fl)
   * 2. Find ObjNums extreme points in a square matrix. Go through selected and Fl, find the fitness value that has the
   * maximum value on the c-th objectives, and fill the c-th coloumn of this matrix with its fitness value. In most
   * cases this matrix is not singular.
   * 3. Compute translated for genes : gene.translatedFitness=gene.fitness-idealPoint
   * 4. Compute translated extreme points : extremePoints.colwise()-=idealPoint
   * 5. The translated extreme points consistis a hyperplane, compute the inverse matrix of extremePoints and find the
   * intercept of this hyperplane. If the matrix is singular, use the diagonal line as intercept.
   * 6. Update the translateFitness : gene.translatedFitness/=intercept.
   */
  void normalize(const std::unordered_set<infoUnit3*>& selected, const std::vector<infoUnitBase_t*>& Fl) const {
    const size_t M = this->objectiveNum();
    stdContainer<const infoUnit3*, ObjNum> extremePtrs;
    heu_initializeSize<ObjNum>::template resize<decltype(extremePtrs)>(&extremePtrs, M);

    Eigen::Array<double, ObjNum, ObjNum> extremePoints;

    Fitness_t ideal, intercepts;

    extremePoints.resize(M, M);
    ideal.resize(M);
    intercepts.resize(M);

    ideal.setConstant(pinfD);

    for (size_t c = 0; c < M; c++) {
      extremePtrs[c] = static_cast<infoUnit3*>(Fl[0]);
      extremePoints.col(c) = extremePtrs[c]->fitnessCache;
    }

    for (auto i : selected) {
      ideal = ideal.min(i->fitnessCache);
      for (size_t objIdx = 0; objIdx < M; objIdx++) {
        if (i->fitnessCache[objIdx] > extremePtrs[objIdx]->fitnessCache[objIdx]) {
          extremePtrs[objIdx] = i;
        }
      }
    }

    for (auto i : Fl) {
      ideal = ideal.min(i->fitnessCache);
      for (size_t objIdx = 0; objIdx < M; objIdx++) {
        if (i->fitnessCache[objIdx] > extremePtrs[objIdx]->fitnessCache[objIdx]) {
          extremePtrs[objIdx] = static_cast<infoUnit3*>(i);
        }
      }
    }

    for (size_t c = 0; c < M; c++) {
      extremePoints.col(c) = extremePtrs[c]->fitnessCache - ideal;
    }

    if (isSingular(extremePoints)) {
      for (size_t r = 0; r < M; r++) {
        intercepts[r] = extremePoints(r, r);
      }
    } else {
      extremePoints2Intercept(extremePoints, &intercepts);
    }

    for (auto i : selected) {
      i->translatedFitness = (i->fitnessCache - ideal) / intercepts;
    }

    for (auto i : Fl) {
      static_cast<infoUnit3*>(i)->translatedFitness = (i->fitnessCache - ideal) / intercepts;
    }
  }  //  normalize

  /**
   * \brief Find the nearest RP of a gene.
   *
   * This word distance doesn't refer to euclidean distance but projection distance from the origin.
   *
   * This step is relatively slow because of great amount of matrix computation.
   *
   * \param s Gene's fitness value
   * \param dist Pass distance through a double pointer.
   * \return size_t The col-index of the nearest RP
   */
  size_t findNearest(const Fitness_t& s, double* dist) const {
    const auto& w = this->referencePoses;
    // w.transpose times s
    auto wT_s = w.matrix().transpose() * s.matrix();
    // w.transpose times s times w
    auto wT_s_w = w.rowwise() * (wT_s.array().transpose());
    // w.transpose times s times w (colwise normalized)
    Eigen::Array<double, ObjNum, Eigen::Dynamic> norm_wTsw = wT_s_w.rowwise() / (w.colwise().squaredNorm());
    auto s_sub_norm_wTsw = norm_wTsw.colwise() - s;
    auto distance = s_sub_norm_wTsw.colwise().squaredNorm();

    int minIdx;
    *dist = distance.minCoeff(&minIdx);
    return minIdx;
  }

  /**
   * \brief Associate a selected ene with a reference point.
   *
   * \param selected hash set of selected genes
   */
  void associate(const std::unordered_set<infoUnit3*>& selected) const {
    for (auto i : selected) {
      i->closestRefPoint = findNearest(i->translatedFitness, &i->distance);
    }
  }

  /**
   * \brief Associate an unselected gene with a reference point.
   *
   * \note Fl_src is a vector of infoUnit3*, while Fl_dst uses index of reference point as its key while gene as
   * value.\n A unordered_multimap is employed cause a RP may be associated by more than one genes, but further
   * procedure requires to find genes by their RP.
   *
   * In this function Fl's associated RP will be found and stored in a unordered_multimap.
   *
   * \param Fl_src Source of Fl
   * \param Fl_dst Destination of Fl
   */
  void associate(const std::vector<infoUnitBase_t*>& Fl_src,
                 std::unordered_multimap<RefPointIdx_t, infoUnit3*>* Fl_dst) const {
    for (auto j : Fl_src) {
      infoUnit3* i = static_cast<infoUnit3*>(j);
      RefPointIdx_t idx = findNearest(i->translatedFitness, &i->distance);
      i->closestRefPoint = idx;
      Fl_dst->emplace(idx, i);
    }
  }  // associate

  /**
   * \brief This niche reservation procedure tends to select genes that are less crowded.
   *
   * In NSGA3, a gene is relatively less crowded means that its RP is associated by fewer genes.
   * 1. The algorithm goes through all reference points, firstly trying to find a RP with least nicheCount (if multiple
   * RPs has the same nicheCount, choose one stochastically).
   * 2. After chooseing a RP, NSGA3 then tries to find its associated genes in Fl(Always remember we need to select part
   * of Fl). If multiple genes are found, choose the closest gene. If no gene is found, erase this RP and find a new
   * one.
   * 3. The gene chosen in previous step is emplaced to selected, and its RP's nicheCount adds by one.
   *
   * Run the previous steps until selected's size is equal to assigend population size.
   *
   * \param selected Selected genes (refers to S_t/F_l in the paper)
   * \param Fl (F_l in the paper)
   * \param refPoints hash map (RP as key and niche count as value)
   */
  void nichePreservation(std::unordered_set<infoUnit3*>* selected,
                         std::unordered_multimap<RefPointIdx_t, infoUnit3*>* Fl,
                         std::unordered_map<RefPointIdx_t, size_t>* refPoints) const {
    for (auto i : *selected) {
      refPoints->operator[](i->closestRefPoint)++;
    }

    std::vector<std::unordered_map<RefPointIdx_t, size_t>::iterator> minNicheIterators;
    minNicheIterators.reserve(refPoints->size());

    std::pair<typename std::unordered_multimap<RefPointIdx_t, infoUnit3*>::iterator,
              typename std::unordered_multimap<RefPointIdx_t, infoUnit3*>::iterator>
        associatedGenesInFl;

    while (selected->size() < this->_option.populationSize) {
      findMinSet(*refPoints, &minNicheIterators);
      auto curRefPoint = minNicheIterators[ei_randIdx(minNicheIterators.size())];
      size_t rhoJ = curRefPoint->second;

      associatedGenesInFl = Fl->equal_range(curRefPoint->first);

      if (associatedGenesInFl.first != associatedGenesInFl.second) {  //  not empty
        typename std::unordered_multimap<RefPointIdx_t, infoUnit3*>::iterator pickedGene;
        if (rhoJ == 0) {
          // find element in associatedGenesInFl with minimum distance
          typename std::unordered_multimap<RefPointIdx_t, infoUnit3*>::iterator minGene = associatedGenesInFl.first;
          for (auto it = associatedGenesInFl.first; it != associatedGenesInFl.second; ++it) {
            if (it->second->distance < minGene->second->distance) {
              minGene = it;
            }
          }
          pickedGene = minGene;
        } else {
          // pick a random member in associatedGenesInFl
          size_t N = 0;
          for (auto it = associatedGenesInFl.first; it != associatedGenesInFl.second; ++it) {
            N++;
          }
          for (auto it = associatedGenesInFl.first; it != associatedGenesInFl.second; ++it) {
            if (ei_randD() * N <= 1) {
              pickedGene = it;
              break;
            }
            N--;
          }
        }

        selected->emplace(pickedGene->second);
        Fl->erase(pickedGene);
        curRefPoint->second++;
      } else {
        refPoints->erase(curRefPoint);
      }
    }  //  end while
  }

  /**
   * \brief Find a subset of Refpoints that has least nicheCount.
   *
   * \param refPoints Refpoints with their nicheCount
   * \param minNicheIterators Destination vector to put the result.
   */
  inline static void findMinSet(std::unordered_map<RefPointIdx_t, size_t>& refPoints,
                                std::vector<std::unordered_map<RefPointIdx_t, size_t>::iterator>* minNicheIterators) {
    minNicheIterators->clear();
    size_t minNiche = 0xFFFFFFFF;
    for (auto i : refPoints) {
      minNiche = std::min(minNiche, i.second);
    }
    for (auto it = refPoints.begin(); it != refPoints.end(); ++it) {
      if (it->second == minNiche) {
        minNicheIterators->emplace_back(it);
      }
    }
  }

 private:
  // make reference points with Das and Dennis’s method recrusively.
  void pri_makeRP(const size_t dimN, const size_t precision, const size_t curDim, const size_t, const size_t accum,
                  Fitness_t* rec, std::vector<Fitness_t>* dst) const {
    if (curDim + 1 >= dimN) {
      rec->operator[](dimN - 1) = 1.0 - double(accum) / precision;
      dst->emplace_back(*rec);
      return;
    }

    for (size_t p = 0; p + accum <= precision; p++) {
      if (curDim >= 0) rec->operator[](curDim) = double(p) / precision;
      pri_makeRP(dimN, precision, curDim + 1, p, accum + p, rec, dst);
    }
  }

  void pri_startRP(const size_t dimN, const size_t precision, std::vector<Fitness_t>* dst) const {
    Fitness_t rec;
#if __cplusplus >= 201703L
    if constexpr (ObjNum == Eigen::Dynamic) {
      rec.resize(this->objectiveNum());
    }
#else
    if (ObjNum == Eigen::Dynamic) {
      rec.resize(this->objectiveNum());
    }
#endif

    pri_makeRP(dimN, precision, 0, 0, 0, &rec, dst);
  }

  inline static bool isSingular(const Eigen::Array<double, ObjNum, ObjNum>& mat) {
    return std::abs(mat.matrix().determinant()) <= 1e-10;
  }

  inline static void extremePoints2Intercept(const Eigen::Array<double, ObjNum, ObjNum>& P, Fitness_t* intercept) {
    auto P_transpose_inv = P.transpose().matrix().inverse();
    auto ONE = Eigen::Matrix<double, ObjNum, 1>::Ones(P.cols(), 1);
    auto one_div_intercept = (P_transpose_inv * ONE).array();
    *intercept = 1.0 / one_div_intercept;
  }
};

#define HEU_MAKE_NSGA3ABSTRACT_TYPES(Base_t)    \
  HEU_MAKE_NSGABASE_TYPES(Base_t)               \
  using RefMat_t = typename Base_t::RefMat_t;   \
  using infoUnit3 = typename Base_t::infoUnit3; \
  using RefPointIdx_t = typename Base_t::RefPointIdx_t;

}  //  namespace internal

}  //  namespace heu

#endif  //  HEU_NSGA3ABSTRACT_HPP
