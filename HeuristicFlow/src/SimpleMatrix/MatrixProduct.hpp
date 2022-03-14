/*
 Copyright © 2022  TokiNoBug
This file is part of Heuristic.

    Heuristic is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Heuristic is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Heuristic.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef Heu_MATRIXPRODUCT_HPP
#define Heu_MATRIXPRODUCT_HPP

#include <type_traits>
#include <assert.h>
#include "../../Global"
#include "MatrixDynamicSize.hpp"
#include "MatrixFixedSize.hpp"
#include "MatrixMap.hpp"

namespace Heu {
/**
 * @brief Matrix of different sizes
 * 
 * @tparam Scalar_t Type of element
 * @tparam rS Row num
 * @tparam cS Coloumn num
 */
template<typename Scalar_t,size_t rS,size_t cS>
using Matrix_t = typename std::conditional<(rS!=Dynamic)&&(cS!=Dynamic),MatrixFixedSize<Scalar_t,rS,cS>,MatrixDynamicSize<Scalar_t>>::type;

template<typename Scalar_t=double,size_t _M=Dynamic,size_t _N=Dynamic,size_t _S=Dynamic>
void MatrixProduct(const Matrix_t<Scalar_t,_M,_N> & A,
                const Matrix_t<Scalar_t,_N,_S> & B,
                Matrix_t<Scalar_t,_M,_S> * dst) {


assert(A.cols()==B.rows());

const size_t M=A.rows();
const size_t N=A.cols();
const size_t S=B.cols();

if
#if __cplusplus >=201703L
    constexpr
#endif
        (_M==Dynamic||_S==Dynamic) {
    dst->resize(M,S);
}

for(size_t i=0;i<M;i++) {
    for(size_t j=0;j<S;j++) {
        Scalar_t sum=0;
        for(size_t k=0;k<N;k++) {
            sum+=A(i,k)*B(k,j);
        }
        dst->operator()(i,j)=sum;
    }
}

}

}//  namespace Heu

   

#endif  //  Heu_MATRIXPRODUCT_HPP
