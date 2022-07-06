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

#ifndef HEU_TESTFUNCTIONS_HPP
#define HEU_TESTFUNCTIONS_HPP

#include "InternalHeaderCheck.h"

#include <type_traits>
#include <HeuristicFlow/Global>

#include "testFunctionsCommon.hpp"

#include "SOFunctions.hpp"
#include "MOFunctions.hpp"

namespace heu {
namespace internal {

struct emptyStruct0 {};
struct emptyStruct1 {};
struct emptyStruct2 {};
struct emptyStruct3 {};
struct emptyStruct4 {};
struct emptyStruct5 {};
struct emptyStruct6 {};
struct emptyStruct7 {};
struct emptyStruct8 {};

template <typename Var_t, class Fitness_t = double, class Arg_t = void>
struct SOFunctions
    : public SOFunctionsX<Var_t, Fitness_t, Arg_t>,
      public std::conditional<sizeMayMatch<Var_t, 2>::value, SOFunctions2<Var_t, Fitness_t, Arg_t>,
                              emptyStruct0>::type {};

template <typename Var_t, class Fitness_t, class Arg_t>
struct MOFunctions
    : public std::conditional<sizeMayMatch<Var_t, 1>::value && sizeMayMatch<Fitness_t, 2>::value,
                              MOFunctions12<Var_t, Fitness_t, Arg_t>, emptyStruct1>::type,
      public std::conditional<sizeMayMatch<Var_t, 2>::value && sizeMayMatch<Fitness_t, 2>::value,
                              MOFunctions22<Var_t, Fitness_t, Arg_t>, emptyStruct2>::type,
      public std::conditional<sizeMayMatch<Var_t, 2>::value && sizeMayMatch<Fitness_t, 3>::value,
                              MOFunctions23<Var_t, Fitness_t, Arg_t>, emptyStruct3>::type,
      public std::conditional<sizeMayMatch<Fitness_t, 2>::value,
                              MOFunctionsX2<Var_t, Fitness_t, Arg_t>, emptyStruct4>::type,
      public std::conditional<sizeMayMatchDTLZ1to7<Var_t, Fitness_t>::value,
                              DTLZ1to7<Var_t, Fitness_t, Arg_t>, emptyStruct5>::type,
      public std::conditional<sizeMayMatchDTLZ89<Var_t, Fitness_t>::value,
                              DTLZ89<Var_t, Fitness_t, Arg_t>, emptyStruct5>::type {
  static_assert((!array_traits<Fitness_t>::isFixedSize) || (array_traits<Fitness_t>::sizeCT >= 2),
                "The size of Fitness_t must be greater than 1.");
};

}  //  namespace internal

template <typename Var_t, class Fitness_t = double, class Arg_t = void>
struct testFunctions
    : public std::conditional<std::is_floating_point<Fitness_t>::value,
                              internal::SOFunctions<Var_t, Fitness_t, Arg_t>,
                              internal::MOFunctions<Var_t, Fitness_t, Arg_t>>::type {};

}  //  namespace heu
#endif
