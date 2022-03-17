// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef Heu_MACROS_HPP
#define Heu_MACROS_HPP

namespace Heu {

///make a member _name in type type
#define Heu_MAKE_MEMBER(type,name) \
    type _##name;


#define Heu_MAKE_VAL_READFUN(type,name) \
    type name() const {return _##name;}

#define Heu_MAKE_VAL_WRITEFUN(type,name) \
    void set_##name(type __##name) {_##name=__##name;}

#define Heu_MAKE_REF_READFUN(type,name) \
    const type & name() const {return _##name;}

#define Heu_MAKE_REF_WRITEFUN(type,name) \
    void set_##name(const type & __##name) {_##name=__##name;}

#define Heu_MAKE_VAL_FUNCTIONS(type,name) \
    Heu_MAKE_VAL_READFUN(type,name) \
    void set_##name(type __##name) {_##name=__##name;}

#define Heu_MAKE_REF_RWFUNCTIONS(type,name) \
    Heu_MAKE_REF_READFUN(type,name) \
    Heu_MAKE_REF_WRITEFUN(type,name)


#define Heu_MAKE_RWATTRIBUTE_REF(type,name) \
protected: \
    Heu_MAKE_MEMBER(type,name) \
public: \
    Heu_MAKE_REF_RWFUNCTIONS(type,name)

#define Heu_MAKE_READONLY_ATTRIBUTE_REF(type,name) \
protected: \
    Heu_MAKE_MEMBER(type,name) \
public: \
    Heu_MAKE_REF_READFUN(type,name)

}

#endif // MACROS_HPP
