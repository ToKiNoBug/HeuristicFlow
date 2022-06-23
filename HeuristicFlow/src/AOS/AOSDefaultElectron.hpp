#ifndef HEU_AOSDEFAULTELECTRON_HPP
#define HEU_AOSDEFAULTELECTRON_HPP

#include "InternalHeaderCheck.hpp"

namespace heu {
namespace internal {
template <typename Var_t, class Fitness_t>
struct DefaultElectron {
  DefaultElectron() { isComputed = false; }
  ~DefaultElectron() = default;
  Var_t state;
  Fitness_t energy;
  bool isComputed;
  inline void setUncomputed() noexcept { isComputed = false; }
};
}  // namespace internal
}  // namespace heu

#endif  //  HEU_AOSDEFAULTELECTRON_HPP