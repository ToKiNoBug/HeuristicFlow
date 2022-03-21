// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_ENUMERATIONS_HPP
#define EIGEN_HEU_ENUMERATIONS_HPP

#include <stdint.h>

#include "InternalHeaderCheck.h"

namespace Eigen {

/**
 * \ingroup HEU_Global
 * \brief Whether to record trainning curve of not
 *
 */
enum RecordOption : uint8_t {
  RECORD_FITNESS = true,       ///< The solver will record fitness values of every generation when running
  DONT_RECORD_FITNESS = false  ///< The solver won't record fitness value.
};

/**
 * \ingroup HEU_Global
 * \brief Convert enumeration to string
 *
 * \param r The enum value
 * \return const char* Name of the value.
 */
inline const char* Enum2String(RecordOption r) {
  switch (r) {
    case RECORD_FITNESS:
      return "RECORD_FITNESS";
    case DONT_RECORD_FITNESS:
      return "DONT_RECORD_FITNESS";
  }
}

/**
 * \ingroup HEU_Global
 * \brief Optimization direction
 *
 */
enum FitnessOption : uint8_t {
  FITNESS_LESS_BETTER = false,    ///< Less fitness value is better
  FITNESS_GREATER_BETTER = true,  ///< Greater fitness value is better
};

/**
 * \ingroup HEU_Global
 * \brief Convert enumeration to string
 *
 * \param f The enum value
 * \return const char* Name of the value.
 */
inline const char* Enum2String(FitnessOption f) {
  switch (f) {
    case FITNESS_LESS_BETTER:
      return "FITNESS_LESS_BETTER";
    case FITNESS_GREATER_BETTER:
      return "FITNESS_GREATER_BETTER";
  }
}

/**
 * \ingroup HEU_Global
 *
 * \brief Which type of container to use.
 *
 * \note Std uses std::vector for dynamic and std::array for fixed.\n\n
 *  While Eigen refers to Eigen::Array<scalar_t,size,1>. If you hope
 *  to use Eigen's matrices or even tensors, inherit it and reload
 *  operator[] to avoid size-checking.\n\n
 *  Custom types should at least be able to act like a vector.
 *  It must has opeartor[](int) that provides random access and
 *  member function size() const that returns the number of elements.
 *  Also if it's dynamic-sized, it should have function resize(int) to change its size.
 */
enum DoubleVectorOption {
  Std = 'S',    ///< C++ standard containers
  Eigen = 'E',  ///< Eigen containers
  Custom = 'C'  ///< User defined custom types.
};

/**
 * \ingroup HEU_Global
 * \brief Convert enumeration to string
 *
 * \param e The enum value
 * \return const char* Name of the value.
 */
inline const char* Enum2String(DoubleVectorOption e) {
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
 * \ingroup HEU_Global
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
 * \ingroup HEU_Global
 * \brief Convert enumeration to string
 *
 * \param b The enum value
 * \return const char* Name of the value.
 */
inline const char* Enum2String(BoxShape b) {
  switch (b) {
    case RECTANGLE_BOX:
      return "Non-square box";
    case SQUARE_BOX:
      return "Square box";
  }
}

/**
 * \ingroup HEU_Global
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
 * \ingroup HEU_Global
 * \brief Convert enumeration to string
 *
 * \param e The enum value
 * \return const char* Name of the value.
 */
inline const char* Enum2String(EncodeType e) {
  switch (e) {
    case Real:
      return "Real encoding";
    case Binary:
      return "Binary encoding";
    case Symbolic:
      return "Symbolic encoding";
  }
}

}  //    namespace Eigen

#endif  // EIGEN_HEU_ENUMERATIONS_HPP
