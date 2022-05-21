

#ifndef HEU_CONSTANTS_HPP
#define HEU_CONSTANTS_HPP

#include <stdint.h>
#include <limits>

#include "InternalHeaderCheck.h"

namespace heu {

namespace internal {

/**
 * \ingroup HEU_GLOBAL
 * \brief Positive infinity float
 *
 */
const float pinfF = std::numeric_limits<float>::infinity();

/**
 * \ingroup HEU_GLOBAL
 * \brief Positive infinity double
 *
 */
const double pinfD = std::numeric_limits<double>::infinity();

/**
 * \ingroup HEU_GLOBAL
 * \brief Negativev infinity float
 *
 */
const float ninfF = -pinfF;

/**
 * \ingroup HEU_GLOBAL
 * \brief Negative infinity double
 *
 */
const double ninfD = -pinfD;

}  //  namespace internal

}  //  namespace heu

#endif  // HEU_CONSTANTS_HPP
