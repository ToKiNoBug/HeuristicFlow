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

#ifndef HEU_ENUMERATIONS_HPP
#define HEU_ENUMERATIONS_HPP

#include <stdint.h>

#include "InternalHeaderCheck.h"

namespace heu {

/**
 * \ingroup HEU_GLOBAL
 * \brief Whether to record trainning curve of not
 *
 * If `RECORD_FITNESS` is assigned to a solver, a specialization of its internal base class will be
 * activated resulting in an extra member `std::vector<Fitness_t> _record`. In that conditon, the
 * solver will record the best fitness of each generation when it's running.
 *
 */
enum RecordOption : uint8_t {
  RECORD_FITNESS =
      true,  ///< The solver will record fitness values of every generation when running
  DONT_RECORD_FITNESS = false  ///< The solver won't record fitness value.
};

/**
 * \ingroup HEU_GLOBAL
 * \brief Convert enumeration to string
 *
 * \param r The enum value
 * \return const char* Name of the value.
 */
inline const char* Enum2String(const RecordOption r) noexcept {
  switch (r) {
    case RECORD_FITNESS:
      return "RECORD_FITNESS";
    case DONT_RECORD_FITNESS:
      return "DONT_RECORD_FITNESS";
  }
}

/**
 * \ingroup HEU_GLOBAL
 * \brief Optimization direction
 *
 * It's vital to tell which direction is better. If `FITNESS_LESS_BETTER` is assigned, a solver will
 * tries the find the minimum value and vise versa.
 */
enum FitnessOption : uint8_t {
  FITNESS_LESS_BETTER = false,    ///< Less fitness value is better
  FITNESS_GREATER_BETTER = true,  ///< Greater fitness value is better
};

/**
 * \ingroup HEU_GLOBAL
 * \brief Convert enumeration to string
 *
 * \param f The enum value
 * \return const char* Name of the value.
 */
inline const char* Enum2String(const FitnessOption f) noexcept {
  switch (f) {
    case FITNESS_LESS_BETTER:
      return "FITNESS_LESS_BETTER";
    case FITNESS_GREATER_BETTER:
      return "FITNESS_GREATER_BETTER";
  }
}

/**
 * \ingroup HEU_GLOBAL
 *
 * \brief Which type of container to use.
 *
 * \note Std uses std::vector for dynamic and std::array for fixed.\n\n
 *  While Eigen refers to `Eigen::Array<scalar_t,size,1>`. If you hope
 *  to use Eigen's matrices or even tensors, inherit it and reload
 *  `operator[](int)` to avoid size-checking.\n\n
 *  Custom types should at least be able to act like a vector.
 *  It must has `opeartor[](int)` that provides random access and
 *  member function `size() const` that returns the number of elements.
 *  Also if it's dynamic-sized, it should have function `resize(int)` to change its size.
 */
enum ContainerOption {
  Std = 'S',    ///< C++ standard containers
  Eigen = 'E',  ///< Eigen containers
  Custom = 'C'  ///< User defined custom types.
};

/**
 * \ingroup HEU_GLOBAL
 * \brief Convert enumeration to string
 *
 * \param e The enum value
 * \return const char* Name of the value.
 */
inline const char* Enum2String(const ContainerOption e) noexcept {
  switch (e) {
    case Std:
      return "C++ std vector/array";
    case Eigen:
      return "Eigen Array";
    default:
      return "Unknown custom types";
  }
}

/**
 * \ingroup HEU_GLOBAL
 * \brief The type of a box-constraint
 *
 * \note Box constraint is the most regular encoding types
 * in evolutionary algoithms.\n
 * Square box is a N-dim hyperbox that each dimension share
 * a same range. And the rest are called Rectangle box or
 * non-square box.
 */
enum BoxShape {
  SQUARE_BOX,    ///< Each dimension has the same range.
  RECTANGLE_BOX  ///< Some dimensions have different ranges with others.
};

/**
 * \ingroup HEU_GLOBAL
 * \brief Convert enumeration to string
 *
 * \param b The enum value
 * \return const char* Name of the value.
 */
inline const char* Enum2String(const BoxShape b) noexcept {
  switch (b) {
    case RECTANGLE_BOX:
      return "Non-square box";
    case SQUARE_BOX:
      return "Square box";
  }
}

/**
 * \ingroup HEU_GLOBAL
 * \brief Encoding type of a box constraint
 *
 * \note Various encoding types are used in evoluntionary
 * algorithms. This enum is usually used together with BoxShape
 *
 */
enum EncodeType {
  Real,    ///< Encode in floating-point numbers.
  Binary,  ///< Encode in binaries, 0 or 1.
  // Integer,

  Symbolic  ///< Encode in symbolic integers.Encode in symbolic integers.
};

/**
 * \ingroup HEU_GLOBAL
 * \brief Convert enumeration to string
 *
 * \param e The enum value
 * \return const char* Name of the value.
 */
inline const char* Enum2String(const EncodeType e) noexcept {
  switch (e) {
    case Real:
      return "Real encoding";
    case Binary:
      return "Binary encoding";
    case Symbolic:
      return "Symbolic encoding";
  }
}

/**
 * \ingroup HEU_GLOBAL
 * \brief Methods of selection in SOGA
 *
 * \note In documents of every method, N refers to the population size BEFORE selection, and K
 * refers to the user-assigned population size (population size AFTER selection). In each select
 * procedure, K genes must be selected from N genes to form the next generation, while the rest N-K
 * genes will be eliminated.
 *
 */
enum SelectMethod : uint32_t {
  /**
   * \brief The most traditional strategy in GA. The probability a gene to be selected is in
   * proportion to its fitness.
   *
   * This method requires fitness values to be of similiar magnitude.
   *
   * \note Usually this strategy requires that the fitness value should always be greater than 0,
   * and a greater fitness value must symbolize a better gene. In the exact implementation, fitness
   * values are firstly inversed if a less fitness value means better(no operation otherwise), and
   * then subtracted with the minimum fitness value.
   *
   *
   */
  RouletteWheel,

  /**
   * \brief The tournament method picks 2 or 3 genes into the "tournament" randomly and select the
   * best one in the "tournament". This process will be repeated many times until the selection
   * finished.
   *
   * This method has no requirement to the fitness value.
   *
   * \note Since all genes in the tournament will be put back into the population after a unit
   * selection, genes may be cpoied for several times.
   *
   * \note This method has a parameter: the tournament size. Its default value is 3. A greater value
   * gives more change of survival for worse genes, which harms the convergence ability but may
   * prevent the solver from being pre-matured.
   */
  Tournament,

  /**
   * \brief The truncation method sort the whole population by fitness value, and the first K genes
   * will be selected, while the rest are eliminated.
   *
   * This method has no requirement to the fitness value.
   *
   * \note In previous versions, SOGA has only one selection method, which is tournament actually.
   * This method usually can't work pretty well even with the simplist testing fucntion.
   *
   */
  Truncation,

  /**
   * \brief The Monte Carlo method select genes purly stochastically. It is neither searching for
   * better nor worse, but ramdonly.
   *
   * This method has no requirement to the fitness value.
   *
   * \note The method is not used for solving problems, but to compare with other methods.
   */
  MonteCarlo,

  /**
   * \brief The probability selection method is a variant of RouletteWheel. In this method, the
   * probability p (refers the the probability of being selected) in RouletteWheel is mapped to
   * p^K*(1-p)^(N-K), and the probability of being selected is in proportion with this mapped
   * probability.
   *
   * This method requires fitness values to be of similiar magnitude.
   *
   */
  Probability,

  /**
   * \brief The LinearRank method sorts the whole population, and then evaluates probability values
   * by ranks of genes. The best and worst genes' probabilities of being selected are assigned by
   * users (usually 0.8 and 0.2 respectively). Probabilities of the rest genes are computed
   * linearly.
   *
   * This method has no requirement to fitness values.
   *
   * \note This method has 2 parameters, which are best and worsts' probabilities of being selected.
   * It must be in range [0,1], and the best must have greater chance of being selected than the
   * worst.
   *
   */
  LinearRank,

  /**
   * \brief The ExponentialRank method works like LinearRank, but probability are evaluated by
   * exponential function. Given a base number c, and the i-th best gene's probability is in
   * proportion with c^idx. If you hope the fitness value to be greater, c should be a positive
   * number, vise versa. c=0 means that all genes have same chance to be selected, no matter how
   * good they are.
   *
   *
   * This method has no requirement to fitness values.
   *
   * \note This method introducted one parameter: the exponential base number c.   *
   */
  ExponentialRank,

  /**
   * \brief The Boltzmann method works by exponential function. The probability of a gene is in
   * proportion with exp(b*f), where b is the selection strength and f is the fitness value.
   *
   * \note This method introduced one parameter: the selection strength b. For maximum problems, b
   * should be positive, while for minimize problems, b should be a negative number.
   */
  Boltzmann,

  /**
   * \brief The SUS(Stochastic Universal Selection) method is a simplification of RouletteWheel. In
   * the RouletteWheel method, the population will be traversed for at least K times, while SUS
   * needs only constant times. It runs faster.
   *
   * This method requires fitness values to be of similiar magnitude.
   */
  StochasticUniversal,

  /**
   * \brief The EliteReserved method applies elitisim to the RouletteWheel method. Several best
   * genes will be selected before the roulette wheel process, which ensures congvergence. This
   * method is proved to be able to converge, which makes it the default choice.
   *
   */
  EliteReserved,

  /**
   * \brief This value indicates that the selection method is determined at runtime. SOGA with this
   * parameter will inherit all 10 selection methods in above, and you can choose which method to
   * use at runtime.
   *
   */
  RunTimeSelectMethod  ///< The selection method is determined at runtime
};

inline const char* Enum2String(const SelectMethod sm) noexcept {
  switch (sm) {
    case RouletteWheel:
      return "RouletteWheel";
    case Tournament:
      return "Tournament";
    case Truncation:
      return "Truncation";
    case MonteCarlo:
      return "MonteCarlo";
    case Probability:
      return "Probability";
    case LinearRank:
      return "LinearRank";
    case ExponentialRank:
      return "ExponentialRank";
    case Boltzmann:
      return "Boltzmann";
    case StochasticUniversal:
      return "StochasticUniversal";
    case EliteReserved:
      return "EliteReserved";
    default:
      return "RunTimeSelectMethod";
  }
}

}  //    namespace heu

#endif  // HEU_ENUMERATIONS_HPP
