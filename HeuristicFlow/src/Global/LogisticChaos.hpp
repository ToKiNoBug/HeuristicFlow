

#ifndef HEU_LOGISTICCHAOS_HPP
#define HEU_LOGISTICCHAOS_HPP

#include <cstdint>
#include <assert.h>

#include "InternalHeaderCheck.h"
#include "Randoms.hpp"

namespace heu {

/**
 * \ingroup HEU_GLOBAL
 * \class LogisticChaos
 * \brief A chaos device to generate logistic chaotic sequence.
 *
 * This device maintance a double and it iterates like this:
 * x=4*x*(1-x)
 *
 */
class LogisticChaos {
 public:
  /**
   * \brief Construct a new Logistic Chaos object with random initial value.
   *
   */
  LogisticChaos() { _value = ei_randD(); }

  /**
   * \brief Construct a new Logistic Chaos object with assigned initial value.
   *
   * \param seed Initial value that should be in range (0,1)
   */
  LogisticChaos(double seed) {
    assert(seed > 0.0);
    assert(seed < 1.0);
    _value = seed;
  }

  /**
   * \brief Iterate and return the next value.
   *
   * \return double Next value.
   */
  double operator()() {
    _value *= miu * (1 - _value);
    return _value;
  }

  /**
   * \brief Get current value.
   *
   * \return double Current value.
   */
  double value() const { return _value; }

  /**
   * \brief Generate a chaotic sequence at given address.
   *
   * \param dst Given address to write chaotic sequence.
   * \param L Lenght of sequence.
   *
   * \note Becareful of memory leaking!
   */
  void makeSequence(double *dst, size_t L) {
    if (L <= 0) return;
    dst[0] = operator()();
    for (size_t i = 1; i < L; i++) {
      dst[i] = miu * dst[i - 1] * (1 - dst[i - 1]);
    }
    _value = dst[L - 1];
  }

  /**
   * \brief Iterate for several times.
   *
   * \param It Times of iterations.
   */
  void iterate(size_t It) {
    for (size_t i = 0; i < It; i++) {
      _value *= miu * (1 - _value);
    }
  }

 private:
  double _value;
  constexpr static const double miu = 4.0;
};

}  //    namespace heu
#endif  // HEU_LOGISTICCHAOS_H
