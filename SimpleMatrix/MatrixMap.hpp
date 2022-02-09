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

#ifndef OptimT_MATRIXMAP_HPP
#define OptimT_MATRIXMAP_HPP

namespace OptimT {

template<typename Scalar_t>
class MatrixMap
{
protected:
using fast_t=typename std::conditional<sizeof(Scalar_t)<=3*sizeof(void*),
    Scalar_t,
    const Scalar_t &>::type;
public:

    MatrixMap(Scalar_t * d,size_t _r,size_t _c) {
        p=d;
        r=_r;
        c=_c;
    }
    /**
     * @brief Construct a new Matrix Map object through shallow copying
     * 
     * @param s 
     */
    MatrixMap(const MatrixMap & s) {
        r=s.r;
        c=s.c;
        p=s.p;
    }

    ~MatrixMap() {
        ;
    }

    void deepCopyTo(MatrixMap * dst) const {
        dst->r=r;
        dst->c=c;
        for(size_t i=0;i<size();i++) {
            dst->operator()(i)=operator()(i);
        }
    }

    inline size_t rows() const {
        return r;
    }

    inline size_t cols() const {
        return c;
    }

    inline size_t size() const {
        return r*c;
    }

    inline Scalar_t * data() const {
        return p;
    }

    inline const Scalar_t * cdata() const {
        return p;
    }

    inline Scalar_t * begin() const {
        return p;
    }

    inline Scalar_t * end() const {
        return p+r*c;
    }

    inline void resize(size_t _r,size_t _c) const {
        r=_r;c=_c;
    }

    inline Scalar_t & operator()(size_t i) const {
        return p[i];
    }

    inline Scalar_t & operator[](size_t i) const {
        return p[i];
    }

    inline Scalar_t & operator()(size_t _r,size_t _c) const {
        return p[r*_c+_r];
    }

    inline void fill(fast_t src) {
        for(auto & i : *this) {
            i=src;
        }
    }

static const bool isFixedSize=false;

protected:
    Scalar_t * p;
    size_t r;
    size_t c;


};


}

#endif  //  OptimT_MATRIXMAP_HPP