/*
 Copyright © 2022  TokiNoBug
This file is part of HeuristicFlow.

    HeuristicFlow is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HeuristicFlow is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HeuristicFlow.  If not, see <https://www.gnu.org/licenses/>.

*/

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
