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

#ifndef HEU_HeuMaths_HPP
#define HEU_HeuMaths_HPP

#include <stdint.h>
#include <type_traits>
#include <cmath>

#include "InternalHeaderCheck.h"

namespace heu {
/*
#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

#ifndef M_PI
#define M_2_PI		0.63661977236758134308
#endif

*/
template <typename A>
inline A square(A a) {
  return a * a;
}

inline double sign(double x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

/**
 * \ingroup HEU_GLOBAL
 * \brief Compute fractorial.
 *
 * \tparam num_t Type of integer
 * \param n An integer no less than 0
 * \return num_t Value of n!
 */
template <typename num_t>
inline num_t fractorial(num_t n) {
  if (n > num_t(1))
    return n * fractorial(n - 1);
  else
    return 1;
}

/**
 * \ingroup HEU_GLOBAL
 * \brief Computer combinatorial number C^K_N.
 * This function is called by NSGA3 when making reference points.
 *
 * \tparam num_t Type of integer
 * \param N To be choosed.
 * \param K Choosed.
 * \return num_t Value of C^K_N.
 */
template <typename num_t>
inline num_t NchooseK(num_t N, num_t K) {
  return fractorial<num_t>(N) / (fractorial<num_t>(K) * fractorial<num_t>(N - K));
}
/*
namespace internal
{

template<typename val_t,int64_t N>
struct Heu_expander

{
static val_t expand(val_t v)
{

return v*Heu_expander<val_t,
std::integral_constant<int64_t,N-((N>0)?(1):(-1))>::value
 >::expand(v);
}
};

template<typename val_t>
struct Heu_expander<val_t,0>

{
static val_t expand(val_t v)
{
return 1;
}
};

} //namespace internal


template<int64_t p,typename val_t>
val_t power(val_t v)
{
if(p>=0)
return internal::Heu_expander<val_t,p>::expand(v);
else
return internal::Heu_expander<val_t,p>::expand(1/v);
}

template<typename val_t>
val_t power(val_t v,int64_t p)
{
double res=1;
if(p<0)
{
p=-p;
v=1/v;
}

for(int64_t i=0;i<p;i++)
{
res*=v;
}
return res;
}

*/

namespace internal {

// Function to implement minimum for multiple inputs
template <typename T>
inline T imp_min(T a, T b) {
  if (a >= b) return b;
  return a;
}

// Function to implement minimum for multiple inputs
template <typename T, typename U, class... Args_t>
inline T imp_min(T a, U b, Args_t... args) {
  static_assert(std::is_same<T, U>::value, "All parameters must be of same types");
  return imp_min(imp_min(a, b), args...);
}

// Function to implement maximum for multiple inputs
template <typename T>
inline T imp_max(T a, T b) {
  if (a <= b) return b;
  return a;
}

// Function to implement maximum for multiple inputs
template <typename T, typename U, class... Args_t>
inline T imp_max(T a, U b, Args_t... args) {
  static_assert(std::is_same<T, U>::value, "All parameters must be of same types");
  return imp_max(imp_max(a, b), args...);
}

}  // namespace internal

/**
 * \ingroup HEU_GLOBAL
 * \brief Find Minimum value for multiple inputs
 *
 * \tparam T Type of input
 * \tparam Args_t They should also be of same type with the first input
 * \param a The first input
 * \param args The rest inputs
 * \return T The minimum value
 */
template <typename T, class... Args_t>
inline T min(T a, Args_t... args) {
  return internal::imp_min(a, args...);
}

/**
 * \ingroup HEU_GLOBAL
 * \brief Find maximum value for multiple inputs
 *
 * \tparam T Type of input
 * \tparam Args_t They should also be of same type with the first input
 * \param a The first input
 * \param args The rest inputs
 * \return T The maximum value
 */
template <typename T, class... Args_t>
inline T max(T a, Args_t... args) {
  return internal::imp_max(a, args...);
}

template <typename scalar_t = double>
inline scalar_t gaussianCurve(const scalar_t x, const scalar_t mu = 0.0,
                              const scalar_t sigma = 1.0) noexcept {
  return 1 / (sigma * std::sqrt(2 * M_PI)) *
         std::exp(-heu::square(x - mu) / (2 * heu::square(sigma)));
}

}  //    namespace heu

#endif  // HEU_Maths_HPP
