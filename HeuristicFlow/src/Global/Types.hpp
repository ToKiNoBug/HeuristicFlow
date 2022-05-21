

#ifndef HEU_TYPES_HPP
#define HEU_TYPES_HPP

#include <Eigen/Core>
#include <type_traits>
#include <array>
#include <vector>

#include "InternalHeaderCheck.h"
#include "Enumerations.hpp"
#include "Constants.hpp"

namespace heu {

/**
 * \ingroup HEU_GLOBAL
 * \brief C++ std container for fixed and dynamic sizes.
 *
 * Use std::vector for dynamic size and std::array as fixed size.
 *
 * \tparam scalar_t Type of element
 * \tparam Dim Size at compile time. Use Eigen::Dynamic for dynamic size.
 */
template <typename scalar_t, int Dim>
using stdContainer =
    typename std::conditional<Dim == Eigen::Dynamic, std::vector<scalar_t>, std::array<scalar_t, Dim> >::type;

/**
 * \ingroup HEU_GLOBAL
 * \brief C++ std container of double.
 *
 * \tparam Size Number of element. Use Eigen::Dynamic for dynamic size.
 */
template <int Size>
using stdVecD_t = stdContainer<double, Size>;

/**
 * \ingroup HEU_GLOBAL
 * \brief Container for different types, sizes and scalar types.
 *
 * \tparam scalar_t Type of elements.
 * \tparam Size Size at compile time. Use Eigen::Dynamic for dynamic size.
 * \tparam DVO Type of container.
 * Use ContainerOption::Eigen for Eigen::Array<scalar_t,Size,1> and
 * other enum values for stdContainer<scalar_t,Size>.
 *
 */
template <typename scalar_t, int Size, ContainerOption DVO>
using Container = typename std::conditional<(DVO != ContainerOption::Eigen), stdContainer<scalar_t, Size>,
                                            Eigen::Array<double, Size, 1> >::type;

// template<ContainerOption dvo,size_t Dim>
// using FitnessVec_t= Container<double,Dim,dvo>;

namespace internal {

template <int _ObjNum>
struct heu_initializeSize {
 public:
  template <typename A>
  inline static void resize(A* v, size_t sz) {
    assert(v->size() == sz);
  }
};

template <>
struct heu_initializeSize<Eigen::Dynamic> {
 public:
  template <typename A>
  inline static void resize(A* v, size_t sz) {
    v->resize(sz);
  }
};

}  //    namespace internal

/**
 * \ingroup HEU_GLOBAL
 * \brief Get the type of element from a type of vector/matrix
 *
 * \tparam Vec_t Type of vector/matrix
 * \return type Type of element
 */
template <class Vec_t>
struct toElement {
  using type = typename std::decay<decltype(Vec_t()[0])>::type;
};

}  //   namespace heu

#endif  // HEU_TYPES_HPP
