// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifndef Heu_MATRIXDYNAMICSIZE_H
#define Heu_MATRIXDYNAMICSIZE_H

#include <stdint.h>
#include <memory>
#include <type_traits>

namespace Heu
{
///col-major matrix with basic types only
template<class Scalar_t,class allocator_t=std::allocator<Scalar_t>>
class MatrixDynamicSize
{
protected:
using fast_t=typename std::conditional<sizeof(Scalar_t)<=3*sizeof(void*),Scalar_t,const Scalar_t &>::type;

using allocT_t = std::allocator_traits<allocator_t>;
public:
    MatrixDynamicSize() {
        rowNum=0;
        colNum=0;
        _capacity=0;
        dataPtr=nullptr;
    };

    MatrixDynamicSize(size_t r,size_t c) {
        rowNum=r;
        colNum=c;
        _capacity=r*c;
        dataPtr=alloc().allocate(_capacity);

            for(size_t i=0;i<_capacity;i++) {
                alloc().construct(dataPtr+i);
            }
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
        if(dataPtr!=nullptr) {
            alloc().deallocate(dataPtr,_capacity);
                for(size_t i=0;i<_capacity;i++) {
                    allocT_t::destroy(alloc(),dataPtr+i);
                }
        }
    };

    using iterator = Scalar_t*;
    using citerator = const Scalar_t*;

    inline iterator begin() noexcept{
        return dataPtr;
    }

    inline iterator end()  noexcept{
        return dataPtr+size();
    }

    inline size_t size() const  noexcept{
        return rowNum*colNum;
    }

    void reserve(size_t s) {
        if(_capacity<=s)
            return;
        
        if(isClass) 
            for(size_t i=0;i<_capacity;i++) {
                alloc().destroy(dataPtr+i);
            }
        alloc().deallocate(dataPtr,_capacity);

        dataPtr=alloc().allocate(s);
        _capacity=s;
        if(isClass) {
            for(size_t i=0;i<_capacity;i++) {
                alloc().construct(dataPtr+i);
            }
        }
    }

    void resize(size_t r,size_t c) {
        if(r*c!=size()) {
            if(r*c>_capacity) {
                if(dataPtr!=nullptr) {
                    if(isClass)
                        for(size_t i=0;i<_capacity;i++) {
                            allocT_t::destroy(alloc(),dataPtr+i);
                        }
                    alloc().deallocate(dataPtr,_capacity);
                }

                dataPtr=alloc().allocate(r*c);
                _capacity=r*c;

                if(isClass)
                    for(size_t i=0;i<_capacity;i++) {
                        allocT_t::construct(alloc(),dataPtr+i);
                    }
            }
            
        }
        
        rowNum=r;
        colNum=c;
    }

    inline size_t rows() const  noexcept{
        return rowNum;
    }

    inline size_t cols() const  noexcept{
        return colNum;
    }

    inline size_t capacity() const  noexcept{
        return _capacity;
    }

    inline const  Scalar_t & operator[](size_t n) const  noexcept{
        return dataPtr[n];
    }

    inline const  Scalar_t & operator()(size_t n) const  noexcept{
        return dataPtr[n];
    }

    inline const  Scalar_t & operator()(size_t r,size_t c) const  noexcept{
        return dataPtr[rowNum*c+r];
    }

    inline Scalar_t & operator[](size_t n) noexcept{
        return dataPtr[n];
    }

    inline Scalar_t & operator()(size_t n) noexcept{
        return dataPtr[n];
    }

    inline Scalar_t & operator()(size_t r,size_t c) noexcept{
        return dataPtr[rowNum*c+r];
    }

    inline Scalar_t * data()  noexcept{
        return dataPtr;
    }

    inline const Scalar_t * data() const  noexcept{
        return dataPtr;
    }

    inline void fill(fast_t src)  noexcept{
        for(auto & i : *this) {
            i=src;
        }
    }

    static const bool isFixedSize=false;

protected:
    Scalar_t * dataPtr;
    size_t rowNum;
    size_t colNum;
    size_t _capacity;

private:
    static const bool isClass=std::is_class<Scalar_t>::value;

    inline static allocator_t & alloc() {
        static allocator_t alloctor;
        return alloctor;
    }
};



}

#endif // MATRIXDYNAMICSIZE_H
