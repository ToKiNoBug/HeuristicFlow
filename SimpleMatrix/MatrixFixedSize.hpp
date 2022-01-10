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

#ifndef MATRIXFIXEDSIZE_H
#define MATRIXFIXEDSIZE_H

#include <stdint.h>

namespace OptimT {

template<class Scalar_t,size_t Rows,size_t Cols>
class MatrixFixedSize
{
public:
    MatrixFixedSize() {};
    ~MatrixFixedSize() {};

    using iterator = Scalar_t*;
    using citerator = const Scalar_t*;

    ///Deep copy
    explicit MatrixFixedSize(const MatrixFixedSize & src) {
        for(size_t i=0;i<size();i++) {
            array[i]=src.array[i];
        }
    }

    iterator begin() {
        return array;
    }

    iterator end() {
        return array+size();
    }

    static size_t size() {
        return Rows*Cols;
    }

    static size_t rows() {
        return Rows;
    }

    static size_t cols() {
        return Cols;
    }

    Scalar_t & operator()(size_t n) {
        return array[n];
    }

    Scalar_t & operator()(size_t r,size_t c) {
        return array[Rows*c+r];
    }

    Scalar_t * data() {
        return array;
    }

    const Scalar_t * cdata() const {
        return array;
    }

    ///Deep copy
    const MatrixFixedSize & operator=(const MatrixFixedSize & src) {
        for(size_t i=0;i<size();i++) {
            array[i]=src.array[i];
        }
        return *this;
    }

protected:
    Scalar_t array[Rows*Cols];
};
}

#endif // MATRIXFIXEDSIZE_H
