// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



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
using Matrix_t = typename std::conditional<(rS!=Runtime)&&(cS!=Runtime),MatrixFixedSize<Scalar_t,rS,cS>,MatrixDynamicSize<Scalar_t>>::type;

template<typename Scalar_t=double,size_t _M=Runtime,size_t _N=Runtime,size_t _S=Runtime>
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
        (_M==Runtime||_S==Runtime) {
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
