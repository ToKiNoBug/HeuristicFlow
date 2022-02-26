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

#ifndef Heu_HeuMaths_HPP
#define Heu_HeuMaths_HPP

#include <stdint.h>
#include <type_traits>
#include <cmath>
namespace Heu {

template<typename A>
inline A square(A a) {
    return a*a;
}

inline double sign(double x) {
    if(x>0) return 1;
    if(x<0) return -1;
    return 0;
}

template<typename num_t>
inline num_t fractorial(num_t n) {
    if(n>num_t(1))
        return n*fractorial(n-1);
    else
        return 1;
}

template<typename num_t>
inline num_t NchooseK(num_t N,num_t K) {
    return fractorial<num_t>(N)/(fractorial<num_t>(K)*fractorial<num_t>(N-K));
}

namespace HeuPrivate {

    template<typename val_t,int64_t N>
    struct expander
    {
        static val_t expand(val_t v) {

            return v*expander<val_t,
                    std::integral_constant<int64_t,N-((N>0)?(1):(-1))>::value
                                                         >::expand(v);
        }
    };

    template<typename val_t>
    struct expander<val_t,0>
    {
        static val_t expand(val_t v) {
            return 1;
        }
    };
    
}   //  namespace HeuPrivate


template<int64_t p,typename val_t>
val_t power(val_t v) {
    if(p>=0)
        return HeuPrivate::expander<val_t,p>::expand(v);
    else
        return HeuPrivate::expander<val_t,p>::expand(1/v);
}

template<typename val_t>
val_t power(val_t v,int64_t p) {
    double res=1;
    if(p<0) {
        p=-p;
        v=1/v;
    }

    for(int64_t i=0;i<p;i++) {
        res*=v;
    }
    return res;
}


namespace HeuPrivate 
{

template<typename T>
inline T imp_min(T a,T b) {
    if (a>=b) {
        return b;
    }
    return a;
}

template<typename T,typename U,class ... Args_t>
inline T imp_min(T a,U b,Args_t ... args) {
    static_assert(std::is_same<T,U>::value,
        "All parameters must be of same types");
    return imp_min(imp_min(a,b),args...);
}

template<typename T>
inline T imp_max(T a,T b) {
    if (a<=b) {
        return b;
    }
    return a;
}

template<typename T,typename U,class ... Args_t>
inline T imp_max(T a,U b,Args_t ... args) {
    static_assert(std::is_same<T,U>::value,
        "All parameters must be of same types");
    return imp_max(imp_max(a,b),args...);
}

}   //  namespace HeuPrivate 
/**
 * @brief minimum value for multiple parameters
 */
template<typename T,class ... Args_t>
inline T min(T a,Args_t ... args) {
    if constexpr (sizeof...(Args_t)>=1)
        return HeuPrivate::imp_min(a,args...);
    else
        return a;
}


/**
 * @brief maximum value for multiple parameters
 */
template<typename T,class ... Args_t>
inline T max(T a,Args_t ... args) {
    if constexpr (sizeof...(Args_t)>=1)
        return HeuPrivate::imp_max(a,args...);
    else
        return a;
}

}

#endif  //  HeuMaths_HPP
