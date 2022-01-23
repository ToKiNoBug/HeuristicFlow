#include <stdint.h>
#include <type_traits>
#include <cmath>
namespace OptimT {

namespace OptimTPrivate {

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
    
}


template<int64_t p,typename val_t>
val_t power(val_t v) {
    if(p>=0)
        return OptimTPrivate::expander<val_t,p>::expand(v);
    else
        return OptimTPrivate::expander<val_t,p>::expand(1/v);
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

}
