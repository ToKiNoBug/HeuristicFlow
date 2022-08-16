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

#ifndef HEU_DEFAULTGENETYPE_HPP
#define HEU_DEFAULTGENETYPE_HPP

#include "InternalHeaderCheck.h"

#include "../../Global"
#include "../../EAGlobal"
#include "GAAbstract.hpp"

#include "IsGene.hpp"

namespace heu {

namespace internal {

/**
 * \ingroup HEU_GENETIC
 * \brief Default type of gene for GA
 *
 * \tparam Var_t Type of decision variable
 * \tparam Fitness_t Type of fitness
 */
template <typename Var_t, typename Fitness_t>
class DefaultGene_t {
 public:
  DefaultGene_t() : is_fitness_computed(false) {
    static_assert(is_GA_gene_v<std::decay_t<decltype(*this)>>);
  }
  Var_t decision_variable;   ///< Value of decision variable
  Fitness_t fitness;         ///< Value of fitness
  bool is_fitness_computed;  ///< Whether the fitness is computed

  ///< Set the fitness to be uncomputed
  inline void set_fitness_uncomputed() noexcept { is_fitness_computed = false; }
};

/**
 * \ingroup HEU_GENETIC
 * \brief Type of gene for MOGAs
 *
 * \tparam Var_t Type of decision variable
 * \tparam N Number of objectives.
 */
template <typename Var_t, int N>
class MOGene_t : public DefaultGene_t<Var_t, Eigen::Array<double, N, 1>> {
 public:
  using Fitness_t = Eigen::Array<double, N, 1>;
};

/**
 * \ingroup HEU_GENETIC
 * \brief Type of gene for NSGA-like algorithm
 *
 * \tparam Var_t Type of decision variable
 * \tparam N Number of objectives.
 */
template <typename Var_t, int N>
class NSGAGene_t : public MOGene_t<Var_t, N> {
 public:
  using Fitness_t = Eigen::Array<double, N, 1>;
  NSGAGene_t() : dominated_by_num(0) {
    static_assert(is_NSGA_gene_v<std::decay_t<decltype(*this)>>);
  }
  /// The number of genes that dominate a gene
  size_t dominated_by_num;
};

/**
 * \ingroup HEU_GENETIC
 * \brief Type of gene for NSGA2.
 *
 * \tparam Var_t Type of decision variable
 * \tparam N Number of objectives.
 */
template <typename Var_t, int N>
class NSGA2Gene_t : public NSGAGene_t<Var_t, N> {
 public:
  NSGA2Gene_t() : congestion(0) { static_assert(is_NSGA2_gene_v<std::decay_t<decltype(*this)>>); }
  double congestion;
};

/**
 * \ingroup HEU_GENETIC
 * \brief Type of gene for NSGA3
 *
 * \tparam Var_t Type of decision variable
 * \tparam N Number of objectives.
 */
template <typename Var_t, int N>
class NSGA3Gene_t : public NSGAGene_t<Var_t, N> {
 public:
  NSGA3Gene_t() : closestRefPoint(0), distance(0) {
    static_assert(is_NSGA3_gene_v<std::decay_t<decltype(*this)>>);
  }
  /// The translated fitness. Normalized fitness is also stored in this.
  typename NSGAGene_t<Var_t, N>::Fitness_t translatedFitness;
  /// The index of its closet RP
  size_t closestRefPoint;
  /// Distance to the closet RP
  double distance;
};

}  // namespace internal

}  // namespace heu

#endif  //  HEU_DEFAULTGENETYPE_HPP