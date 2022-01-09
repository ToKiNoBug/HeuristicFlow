#ifndef SIMPLEMATRIX_H
#define SIMPLEMATRIX_H

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

    Scalar_t & operator()(size_t n) {
        return dataPtr[n];
    }

    Scalar_t & operator()(size_t r,size_t c) {
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

#endif // SIMPLEMATRIX_H
