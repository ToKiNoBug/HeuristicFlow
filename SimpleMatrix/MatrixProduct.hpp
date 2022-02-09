#ifndef OptimT_MATRIXPRODUCT_HPP
#define OptimT_MATRIXPRODUCT_HPP

#include "../OptimTemplates/Global"
#include <type_traits>
#include <assert.h>
#include "MatrixDynamicSize.hpp"
#include "MatrixFixedSize.hpp"
#include "MatrixMap.hpp"

namespace OptimT {
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

if constexpr(_M==Dynamic||_S==Dynamic) {
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

}//  namespace OptimT

   

#endif  //  OptimT_MATRIXPRODUCT_HPP