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

#ifndef HEU_MACROS_HPP
#define HEU_MACROS_HPP

#include <assert.h>
#include <cstdlib>

#define HEU_RELOAD_MEMBERFUCTION_RUN                                         \
  inline void run() noexcept {                                               \
    this->template __impl_run<typename std::decay<decltype(*this)>::type>(); \
  }

#define HEU_DISPLINE \
  ::std::cout << "File : " << __FILE__ << " , Line : " << __LINE__ << ::std::endl;

#ifdef __cplusplus
#define HEU_CPP_STANDARD __cplusplus
#else
#define HEU_CPP_STANDARD _MSC_VER
#endif

#define HEU_ASSERT(expression) \
  {                            \
    assert(expression);        \
    if (!(expression)) {       \
      abort();                 \
    }                          \
  }

#endif  //  HEU_MACROS_HPP