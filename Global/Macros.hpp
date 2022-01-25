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

#ifndef MACROS_HPP
#define MACROS_HPP

namespace OptimT {

///make a member _name in type type
#define OptimT_MAKE_MEMBER(type,name) \
    type _##name;


#define OptimT_MAKE_VAL_READFUN(type,name) \
    type name() const {return _##name;}

#define OptimT_MAKE_VAL_WRITEFUN(type,name) \
    void set_##name(type __##name) {_##name=__##name;}

#define OptimT_MAKE_REF_READFUN(type,name) \
    const type & name() const {return _##name;}

#define OptimT_MAKE_REF_WRITEFUN(type,name) \
    void set_##name(const type & __##name) {_##name=__##name;}

#define OptimT_MAKE_VAL_FUNCTIONS(type,name) \
    OptimT_MAKE_VAL_READFUN(type,name) \
    void set_##name(type __##name) {_##name=__##name;}

#define OptimT_MAKE_REF_RWFUNCTIONS(type,name) \
    OptimT_MAKE_REF_READFUN(type,name) \
    OptimT_MAKE_REF_WRITEFUN(type,name)


#define OptimT_MAKE_RWATTRIBUTE_REF(type,name) \
protected: \
    OptimT_MAKE_MEMBER(type,name) \
public: \
    OptimT_MAKE_REF_RWFUNCTIONS(type,name)

#define OptimT_MAKE_READONLY_ATTRIBUTE_REF(type,name) \
protected: \
    OptimT_MAKE_MEMBER(type,name) \
public: \
    OptimT_MAKE_REF_READFUN(type,name)

}

#endif // MACROS_HPP
