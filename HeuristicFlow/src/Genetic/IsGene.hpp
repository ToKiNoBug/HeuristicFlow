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

#ifndef HEU_ISGENE_HPP
#define HEU_ISGENE_HPP

#include "../../Global"

namespace heu {

namespace {
template <class G>
static inline decltype(G().decision_variable, G().fitness, G().is_fitness_computed,
                       G().set_fitness_uncomputed(), int())
__impl_fun_is_gene(const G &) {
  return 0;
}

static inline void __impl_fun_is_gene(...) { return; }
}  // namespace

template <class G>
constexpr bool is_GA_gene_v = std::is_same_v<decltype(__impl_fun_is_gene(G())), int>;

namespace {
template <class G>
static inline decltype(std::enable_if_t<is_GA_gene_v<G>>(), G().dominated_by_num, int())
__impl_fun_is_NSGA_gene(const G &) {
  return 0;
}

static inline void __impl_fun_is_NSGA_gene(...) { return; }
}  // namespace

template <class G>
constexpr bool is_NSGA_gene_v = std::is_same_v<decltype(__impl_fun_is_NSGA_gene(G())), int>;

namespace {
template <class G>
static inline decltype(std::enable_if_t<is_NSGA_gene_v<G>>(), G().congestion, int())
__impl_fun_is_NSGA2_gene(const G &) {
  return 0;
}

static inline void __impl_fun_is_NSGA2_gene(...) { return; }
}  // namespace

template <class G>
constexpr bool is_NSGA2_gene_v = std::is_same_v<decltype(__impl_fun_is_NSGA2_gene(G())), int>;

namespace {
template <class G>
static inline decltype(std::enable_if_t<is_NSGA_gene_v<G>>(), G().translatedFitness,
                       G().closestRefPoint, G().distance, int())
__impl_fun_is_NSGA3_gene(const G &) {
  return 0;
}

static inline void __impl_fun_is_NSGA3_gene(...) { return; }
}  // namespace

template <class G>
constexpr bool is_NSGA3_gene_v = std::is_same_v<decltype(__impl_fun_is_NSGA3_gene(G())), int>;

}  // namespace heu

#endif  // HEU_ISGENE_HPP