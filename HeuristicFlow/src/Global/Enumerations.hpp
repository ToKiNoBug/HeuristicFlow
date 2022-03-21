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
enum RecordOption : uint8_t { RECORD_FITNESS = true, DONT_RECORD_FITNESS = false };

/// convert enumeration to string
inline const char* Enum2String(RecordOption r) {
  switch (r)

  {
    case RECORD_FITNESS:
      return "RECORD_FITNESS";
    case DONT_RECORD_FITNESS:
      return "DONT_RECORD_FITNESS";
  }
}

/// optimization direction
enum FitnessOption : uint8_t {
  FITNESS_LESS_BETTER = false,
  FITNESS_GREATER_BETTER = true,
};

/// convert enumeration to string
inline const char* Enum2String(FitnessOption f) {
  switch (f)

  {
    case FITNESS_LESS_BETTER:
      return "FITNESS_LESS_BETTER";
    case FITNESS_GREATER_BETTER:
      return "FITNESS_GREATER_BETTER";
  }
}

/// which type of vector to use
enum DoubleVectorOption { Std = 'S', Eigen = 'E', Custom = 'C' };
/// convert enumeration to string
inline const char* Enum2String(DoubleVectorOption e) {
  switch (e)

  {
    case Std:
      return "C++ std vector/array";
    case Eigen:
      return "Eigen Array";
    default:
      return "Unknown custom types";
  }
}

enum BoxShape { RECTANGLE_BOX, SQUARE_BOX };

inline const char* Enum2String(BoxShape b) {
  switch (b) {
    case RECTANGLE_BOX:
      return "Non-square box";
    case SQUARE_BOX:
      return "Square box";
  }
}

enum EncodeType {
  Real,
  Binary,
  // Integer,
  Symbolic
};

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
