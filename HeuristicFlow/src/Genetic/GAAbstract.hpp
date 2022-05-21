
#ifndef HEU_GAABSTRACT_HPP
#define HEU_GAABSTRACT_HPP

#include <type_traits>

#include "InternalHeaderCheck.h"
#include <HeuristicFlow/Global>

namespace heu {

namespace internal {

HEU_MAKE_FUNAREA(iFun, GA)
HEU_MAKE_FUNAREA(mFun, GA)
HEU_MAKE_FUNAREA(fFun, GA)
HEU_MAKE_FUNAREA(cFun, GA)

/**
 * \ingroup HEU_GENETIC
 * \class GAAbstract
 * \brief Internal base class.
 *
 * This class defines function pointer types of the four operator(initialization, computing fitness, crossover and
 * mutation)
 *
 * \tparam Var_t Type of decision variable
 * \tparam Fitness_t Type of fitness
 * \tparam Args_t Type of global parameters
 */
template <typename Var_t, typename Fitness_t, class Args_t>
class GAAbstract {
 public:
  ~GAAbstract(){};
  /// Function to initialize Var
  using initializeFun = void (*)(Var_t *, const Args_t *);
  /// Function to calculate fitness for Var
  using fitnessFun = void (*)(const Var_t *, const Args_t *, Fitness_t *);
  /// Function to apply crossover for Var
  using crossoverFun = void (*)(const Var_t *, const Var_t *, Var_t *, Var_t *, const Args_t *);
  /// Function to apply mutate for Var
  using mutateFun = void (*)(const Var_t *, Var_t *, const Args_t *);
  /// Type alias of other parameters
  using ArgsType = Args_t;

  /**
   * \class iFunBody
   * \brief Base class that maintains the initialization function ptr
   *
   * \tparam i Initialization function ptr.
   */
  template <initializeFun i>
  using iFunBody = typename iFunArea_GA<Var_t *, const Args_t *>::template funBody<i>;

  /**
   * \class fFunBody
   * \brief Base class that maintains the fitness function ptr
   *
   * \tparam f Fitness function ptr.
   */
  template <fitnessFun f>
  using fFunBody = typename fFunArea_GA<const Var_t *, const Args_t *, Fitness_t *>::template funBody<f>;

  /**
   * \class cFunBody
   * \brief Base class that maintains the crossover function ptr
   *
   * \tparam c Crossover function ptr.
   */
  template <crossoverFun c>
  using cFunBody =
      typename cFunArea_GA<const Var_t *, const Var_t *, Var_t *, Var_t *, const Args_t *>::template funBody<c>;

  /**
   * \class mFunBody
   * \brief Base class that maintains the mutation function ptr
   *
   * \tparam m Mutation function ptr.
   */
  template <mutateFun m>
  using mFunBody = typename mFunArea_GA<const Var_t *, Var_t *, const Args_t *>::template funBody<m>;

  const Args_t &args() const { return _args; }  ///< Const reference to the other parameters

  void setArgs(const Args_t &a) { _args = a; }  ///< Set the value of other parameters

  static const bool HasParameters = true;  ///< This member denotes that this Args_t is not void

 protected:
  Args_t _args;  ///< Pesudo global variable that stored inside the solver
};

// Partial specialization when Args_t is void
template <typename Var_t, typename Fitness_t>
class GAAbstract<Var_t, Fitness_t, void> {
 public:
  ~GAAbstract(){};

  /// Function to initialize Var
  using initializeFun = void (*)(Var_t *);
  /// Function to calculate fitness for Var
  using fitnessFun = void (*)(const Var_t *, Fitness_t *);
  /// Function to apply crossover for Var
  using crossoverFun = void (*)(const Var_t *, const Var_t *, Var_t *, Var_t *);
  /// Function to apply mutate for Var
  using mutateFun = void (*)(const Var_t *, Var_t *);

  using ArgsType = void;

  template <initializeFun i>
  using iFunBody = typename iFunArea_GA<Var_t *>::template funBody<i>;

  template <fitnessFun f>
  using fFunBody = typename fFunArea_GA<const Var_t *, Fitness_t *>::template funBody<f>;

  template <crossoverFun c>
  using cFunBody = typename cFunArea_GA<const Var_t *, const Var_t *, Var_t *, Var_t *>::template funBody<c>;

  template <mutateFun m>
  using mFunBody = typename mFunArea_GA<const Var_t *, Var_t *>::template funBody<m>;

  static const bool HasParameters = false;
};

#define HEU_MAKE_GAABSTRACT_TYPES(Base_t)               \
  using initializeFun = typename Base_t::initializeFun; \
  using fitnessFun = typename Base_t::fitnessFun;       \
  using crossoverFun = typename Base_t::crossoverFun;   \
  using mutateFun = typename Base_t::mutateFun;         \
  using ArgsType = typename Base_t::ArgsType;

}  //  namespace internal

}  //  namespace heu

#endif  //  HEU_GAABSTRACT_HPP
