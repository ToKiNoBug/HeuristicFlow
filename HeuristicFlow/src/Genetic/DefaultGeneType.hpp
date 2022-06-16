#ifndef HEU_DEFAULTGENETYPE_HPP
#define HEU_DEFAULTGENETYPE_HPP

#include "InternalHeaderCheck.h"

#include "../../Global"
#include "../../EAGlobal"
#include "GAAbstract.hpp"

namespace heu {

namespace internal {
/**
 * \ingroup HEU_GENETIC
 * \brief Default type of gene for GA
 *
 * \tparam Var_t
 * \tparam Fitness_t
 */
template <typename Var_t, typename Fitness_t>
class DefaultGene_t {
 public:
  using fastFitness_t =
      typename std::conditional<(sizeof(Fitness_t) > sizeof(double)), const Fitness_t &, Fitness_t>::type;

  Var_t self;          ///< Value of decision variable
  Fitness_t _Fitness;  ///< Value of fitness
  bool _isCalculated;  ///< Whether the fitness is computed

  inline bool isCalculated() const noexcept { return _isCalculated; }  ///< If the fitness is computed
  inline void setUncalculated() noexcept { _isCalculated = false; }    ///< Set the fitness to be uncomputed
  inline fastFitness_t fitness() const noexcept { return _Fitness; }   ///< Get fitness
};

}  // namespace internal

}  // namespace heu

#endif  //  HEU_DEFAULTGENETYPE_HPP