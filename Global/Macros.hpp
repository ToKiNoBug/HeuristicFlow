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
