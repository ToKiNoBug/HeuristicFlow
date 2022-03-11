#ifndef FUNBODY_HPP
#define FUNBODY_HPP

#include <type_traits>

namespace Heu
{
#define Heu_MAKE_FUNAREA(funFlag,FunFlag,Suffix) \
template<typename return_t=void,typename ... a> \
class funFlag##Area_##Suffix \
{ \
public: \
    using funPtr_t = return_t(*)(a...); \
private: \
    template<return_t(*_fun)(a...)> \
    class funBodyCT \
    { \
    public: \
        inline constexpr funPtr_t funFlag() const { \
            return _fun; \
        } \
        inline return_t run##FunFlag(a... _a) const { \
            if constexpr (std::is_same<return_t,void>::value) \
                    _fun(_a...); \
            else \
                return _fun(_a...); \
        } \
    private: \
        static_assert (_fun!=nullptr, "Template function mustn't be nullptr"); \
    }; \
     \
    class funBodyRT \
    { \
    public: \
        funBodyRT() { \
            _funPtr=nullptr; \
        } \
        inline funPtr_t funFlag() const { \
            return _funPtr; \
        } \
        inline return_t run##FunFlag(a... _a) const { \
            if constexpr (std::is_same<return_t,void>::value) \
                    _funPtr(_a...); \
            else \
                return _funPtr(_a...); \
        } \
        inline void set##FunFlag(funPtr_t __) { \
            _funPtr=__; \
        } \
    protected: \
        funPtr_t _funPtr; \
    }; \
public: \
    template<return_t(*_fun)(a...)> \
    using funBody= \
        typename std::conditional<_fun==nullptr,funBodyRT,funBodyCT<_fun>>::type; \
};

}   //  namespace Heu

#endif // FUNBODY_HPP
