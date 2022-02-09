/*
 Copyright © 2022  TokiNoBug
This file is part of OptimTemplates.

    OptimTemplates is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OptimTemplates is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OptimTemplates.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef OptimT_STUPIDLU_HPP
#define OptimT_STUPIDLU_HPP
#include "../OptimTemplates/Global"
#include <type_traits>
#include <assert.h>
#include "MatrixDynamicSize.hpp"
#include "MatrixFixedSize.hpp"
#include "MatrixMap.hpp"
namespace OptimT {

/**
 * @brief A simple LU implementation. 
 *      This function is written only to enable NSGA3 to run without Eigen. 
 *      For high-performance inverse matrix, please use Eigen.
 * 
 * @tparam Mat_t Type of matrix
 * @param src source of matrix
 * @param dst result
 */

template<typename Scalar_t,size_t N>
using SquareMat_t = typename std::conditional<N==Dynamic,MatrixDynamicSize<Scalar_t>,MatrixFixedSize<Scalar_t,N,N>>::type;

template<typename Scalar_t,size_t Size>
void SquareMatrixProduct(const SquareMat_t<Scalar_t,Size> & A,
    const SquareMat_t<Scalar_t,Size> & B,
    SquareMat_t<Scalar_t,Size> * res) {

const size_t N=A.rows();
if constexpr (Size==Dynamic) {
res->resize(N,N);
}

for(size_t r=0;r<N;r++) {
    for(size_t c=0;c<N;c++) {
        Scalar_t sum=0;
        for(size_t k=0;k<N;k++) {
            sum+=A(k,c)*B(r,k);
        }
        res->operator()(r,c)=sum;
    }
}

}

template<typename Scalar_t,size_t Size>
void InverseMatrix_LU(const SquareMat_t<Scalar_t,Size> & A,
                      SquareMat_t<Scalar_t,Size> * invA) {

assert(A.rows()==A.cols());


SquareMat_t<Scalar_t,Size> L,U,iL,iU;
const size_t N=A.rows();
if constexpr (Size==Dynamic) {
L.resize(N,N);
U.resize(N,N);
iL.resize(N,N);
iU.resize(N,N);
}


L.fill(0);
U.fill(0);


for(size_t i=0;i<N;i++) {
    U(0,i)=A(0,i);
    L(i,0)=A(i,0)/U(0,0);
    L(i,i)=1;
}


///compute L and U
for(size_t r=0;r<N;r++) {
    for(size_t j=r;j<N;j++) {
        Scalar_t temp=0;
        for(size_t k=0;k<r;k++) {
            temp+=L(r,k)*U(k,j);
        }
        U(r,j)=A(r,j)-temp;
    }

    for(size_t i=r+1;i<N;i++) {
        Scalar_t  temp=0;
        for(size_t k=0;k<r;k++) {
            temp+=L(i,k)*U(k,r);
        }
        L(i,r)=(A(i,r)-temp)/U(r,r);
    }

}

iL=L;
iU=U;
///compute iL and iU
for(size_t j=0;j<N;j++) {
    for(size_t i=j-1;i!=-1;i--) {
        Scalar_t temp=0;
        for(size_t k=i+1;k<=j;k++) {
            temp+=iL(j,k)*L(k,i);
        }
        iL(j,i)=-iL(i,i)*temp;
    }
    iL(j,j)=1/L(j,j);
}

for(size_t j=N-1;j!=-1;j--) {
    iU(j,j)=1/U(j,j);
    for(size_t i=j+1;i<N;i++) {
        Scalar_t temp=0;
        for(size_t k=j+1;k<=i;k++) {
            temp+=U(j,k)*iU(k,i);
        }
        iU(j,i)=-temp/U(j,j);
    }
}

SquareMatrixProduct<Scalar_t,Size>(iL,iU,invA);

}

}   //  namespace OptimT

#endif  //  OptimT_STUPIDLU_HPP
