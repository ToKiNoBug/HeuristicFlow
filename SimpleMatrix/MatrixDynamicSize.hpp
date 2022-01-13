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

#ifndef MATRIXDYNAMICSIZE_H
#define MATRIXDYNAMICSIZE_H

#include <stdint.h>

namespace OptimT
{
///col-major matrix with basic types only
template<class Scalar_t>
class MatrixDynamicSize
{
public:
    MatrixDynamicSize() {
        rowNum=0;
        colNum=0;
        dataPtr=nullptr;
    };

    MatrixDynamicSize(size_t r,size_t c) {
        rowNum=r;
        colNum=c;
        dataPtr=new Scalar_t[r*c];
    }

    ///Deep copy function
    explicit MatrixDynamicSize(const MatrixDynamicSize & src) {
        resize(src.rowNum,src.colNum);
        for(size_t i=0;i<size();i++) {
            dataPtr[i]=src.dataPtr[i];
        }
    }

    ///Deep copy
    const MatrixDynamicSize & operator=(const MatrixDynamicSize & src) {
        resize(src.rowNum,src.colNum);
        for(size_t i=0;i<size();i++) {
            dataPtr[i]=src.dataPtr[i];
        }
        return *this;
    }

    ~MatrixDynamicSize() {
        if(dataPtr!=nullptr)
            delete [] dataPtr;
    };

    using iterator = Scalar_t*;
    using citerator = const Scalar_t*;

    iterator begin() {
        return dataPtr;
    }

    iterator end() {
        return dataPtr+size();
    }

    size_t size() const {
        return rowNum*colNum;
    }

    void resize(size_t r,size_t c) {
        if(r*c!=size())
        {
            if(dataPtr!=nullptr)
                delete [] dataPtr;
            rowNum=r;
            colNum=c;
            dataPtr=new Scalar_t[r*c];
        }
        else {
            //conservative resize
            rowNum=r;
            colNum=c;
        }
    }

    size_t rows() const {
        return rowNum;
    }

    size_t cols() const {
        return colNum;
    }

    Scalar_t & operator()(size_t n) const {
        return dataPtr[n];
    }

    Scalar_t & operator()(size_t r,size_t c) const {
        return dataPtr[rowNum*c+r];
    }

    Scalar_t * data() {
        return dataPtr;
    }

    const Scalar_t * cdata() const {
        return dataPtr;
    }

protected:
    Scalar_t * dataPtr;
    size_t rowNum;
    size_t colNum;



};



}

#endif // MATRIXDYNAMICSIZE_H
