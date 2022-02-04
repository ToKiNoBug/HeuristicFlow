#ifndef OptimT_GAABSTRACT_HPP
#define OptimT_GAABSTRACT_HPP

namespace OptimT {

template<typename Var_t,typename Fitness_t,class Args_t>
class GAAbstract
{
public:
   GAAbstract() {};
    virtual ~GAAbstract() {};

   ///Function to initialize Var
   using initializeFun = void(*)(Var_t*,const Args_t*);
   ///Function to calculate fitness for Var
   using fitnessFun = void(*)(const Var_t*,const Args_t*,Fitness_t*);
   ///Function to apply crossover for Var
   using crossoverFun = void(*)(const Var_t*,const Var_t*,Var_t*,Var_t*,const Args_t*);
   ///Function to apply mutate for Var
   using mutateFun = initializeFun;

    const Args_t & args() const {
        return _args;
    }

    void setArgs(const Args_t & a) {
        _args=a;
    }

    inline void doInitialize(initializeFun iFun,Var_t * v) {
        iFun(v,&_args);
    }

    inline void doFitness(fitnessFun fFun,const Var_t*v,Fitness_t*f) {
        fFun(v,&_args,f);
    }

    inline void doCrossover(crossoverFun cFun,
                                   const Var_t*p1,const Var_t*p2,Var_t*c1,Var_t*c2) {
        cFun(p1,p2,c1,c2,&_args);
    }
    inline void doMutate(mutateFun mFun,Var_t * v) {
        mFun(v,&_args);
    }

    static const bool HasParameters=true;

protected:
Args_t _args;
};


template<typename Var_t,typename Fitness_t>
class GAAbstract<Var_t,Fitness_t,void>
{
public:
    GAAbstract() {};
    virtual ~GAAbstract() {};

    ///Function to initialize Var
    using initializeFun = void(*)(Var_t*);
    ///Function to calculate fitness for Var
    using fitnessFun = void(*)(const Var_t*,Fitness_t*);
    ///Function to apply crossover for Var
    using crossoverFun = void(*)(const Var_t*,const Var_t*,Var_t*,Var_t*);
    ///Function to apply mutate for Var
    using mutateFun = initializeFun;

    inline static void doInitialize(initializeFun iFun,Var_t * v) {
        iFun(v);
    }

    inline static void doFitness(fitnessFun fFun,const Var_t*v,Fitness_t*f) {
        fFun(v,f);
    }

    inline static void doCrossover(crossoverFun cFun,
                                   const Var_t*p1,const Var_t*p2,Var_t*c1,Var_t*c2) {
        cFun(p1,p2,c1,c2);
    }
    inline static void doMutate(mutateFun mFun,Var_t * v) {
        mFun(v);
    }


    static const bool HasParameters=false;
};

#define OptimT_MAKE_GAABSTRACT_TYPES \
using initializeFun = typename Base_t::initializeFun; \
using fitnessFun = typename Base_t::fitnessFun; \
using crossoverFun = typename Base_t::crossoverFun; \
using mutateFun = typename Base_t::mutateFun;

}   //  OptimT

#endif  //  OptimT_GAABSTRACT_HPP
