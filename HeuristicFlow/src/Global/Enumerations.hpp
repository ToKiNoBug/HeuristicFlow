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
inline const char* Enum2String(RecordOption r) noexcept {
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
inline const char* Enum2String(FitnessOption f) noexcept {
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
inline const char* Enum2String(ContainerOption e) noexcept {
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
inline const char* Enum2String(BoxShape b) noexcept {
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
inline const char* Enum2String(EncodeType e) noexcept {
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
 * \brief Methods of selection in GA
 *
 * \note Many selection methods have requirements to the fitness value.
 */
enum SelectMethod {
  RouletteWheel,        ///< Requirements: fitness values are of similiar magnitude
  Tournament,           ///< Requirements: No
  Truncation,           ///< Requirements: No
  MonteCarlo,           ///< Requirements: No
  Probability,          ///< Requirements:
  LinearRank,           ///< Requirements:
  ExponentialRank,      ///< Requirements:
  Boltzmann,            ///< Requirements:
  StochasticUniversal,  ///< Requirements:
  EliteReserved,        ///< Elitism + the RouletteWheel method. Requirements:
  RunTimeSelectMethod   ///< The selection method is determined at runtime
};

}  //    namespace heu

#endif  // HEU_ENUMERATIONS_HPP
