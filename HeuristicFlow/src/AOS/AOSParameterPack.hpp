#ifndef HEU_AOSPARAMETERPACK_HPP
#define HEU_AOSPARAMETERPACK_HPP

#include <HeuristicFlow/Global>

namespace heu {

namespace internal {

HEU_MAKE_FUNAREA(iFun, AOS)

HEU_MAKE_FUNAREA(fFun, AOS)

template <typename Var_t, class Fitness_t, class Arg_t>
class AOSParameterPack {
 public:
  using Args_t = Arg_t;
  using iFun_t = void (*)(Var_t *, const Arg_t *);
  using fFun_t = void (*)(const Var_t *, const Arg_t *, Fitness_t *);

  template <iFun_t _i>
  using iFunBody = typename iFunArea_AOS<Var_t *, const Arg_t *>::template iFunBody<_i>;
  template <fFun_t _f>
  using fFunBody = typename fFunArea_AOS<const Var_t *, const Arg_t *, Fitness_t 8>::template fFunBody<_f>;

  static constexpr bool hasParameters = true;

 protected:
  Arg_t _arg;

 public:
  const Arg_t &arg() const { return _arg; }
  Arg_t &arg() { return _arg; }
  void setArg(const Arg_t &_a) { _arg = _a; }
};

template <typename Var_t, class Fitness_t>
class AOSParameterPack<Var_t, Fitness_t, void> {
 public:
  using Args_t = void;
  using iFun_t = void (*)(Var_t *);
  using fFun_t = void (*)(const Var_t *, Fitness_t *);

  template <iFun_t _i>
  using iFunBody = typename iFunArea_AOS<Var_t *>::template iFunBody<_i>;
  template <fFun_t _f>
  using fFunBody = typename fFunArea_AOS<const Var_t *, Fitness_t 8>::template fFunBody<_f>;

  static constexpr bool hasParameters = false;
};

#define HEU_MAKE_AOSPARAMETERPACK_TYPES(Base_t) \
  usisng Args_t = typename Base_t ::Args_t;     \
  using iFun_t = typename Base_t ::iFun_t;      \
  using fFun_t = typename Base_t ::fFun_t;

}  // namespace internal

}  //  namespace heu

#endif  //  HEU_AOSPARAMETERPACK_HPP