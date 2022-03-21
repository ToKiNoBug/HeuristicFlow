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

/// whether to record trainning curve of not

/**
 * \brief Whether to record trainning curve of not
 *
 */
enum RecordOption : uint8_t {
  /**
   * \brief The solver will record fitness values of every generation when running
   *
   */
  RECORD_FITNESS = true,
  /**
   * \brief The solver won't record fitness value.
   *
   */
  DONT_RECORD_FITNESS = false
};

/**
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
 * \brief Optimization direction
 *
 */
enum FitnessOption : uint8_t {
  /**
   * \brief Less fitness value is better
   *
   */
  FITNESS_LESS_BETTER = false,
  /**
   * \brief Greater fitness value is better
   *
   */
  FITNESS_GREATER_BETTER = true,
};

/**
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

/// which type of vector to use

/**
 * \brief Which type of container to use.
 *
 */
enum DoubleVectorOption {
  /**
   * \brief C++ standard containers (std::vector for dynamic and std::array for fixed)
   *
   */
  Std = 'S',
  /**
   * \brief Eigen containers (Eigen::Array<scalar_t,size,1>)
   *
   */
  Eigen = 'E',

  /**
   * \brief Use's custom types
   *
   */
  Custom = 'C'
};

/**
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
 * \brief The type of a box-constraint
 *
 */
enum BoxShape {
  /**
   * \brief The box is a square box, which means it has the same range in every dimensions.
   *
   */
  SQUARE_BOX,
  /**
   * \brief A non-square box don't a same range in every dimensions.
   *
   */
  RECTANGLE_BOX
};

/**
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
 * \brief Encoding type of a box constraint
 *
 */
enum EncodeType {
  /**
   * \brief Encode in floating-point numbers.
   *
   */
  Real,
  /**
   * \brief Encode in binaries, 0 or 1.
   *
   */
  Binary,
  // Integer,

  /**
   * \brief Encode in symbolic integers.
   *
   */
  Symbolic
};

/**
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
