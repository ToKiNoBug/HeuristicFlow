
#ifndef HEU_MOGABASE_HPP
#define HEU_MOGABASE_HPP

#include <assert.h>

#include "InternalHeaderCheck.h"
#include "MOGAAbstract.hpp"

#ifndef HEU_MAX_RUNTIME_OBJNUM

/**
 * \ingroup HEU_GENETIC
 * \brief Maximum runtime objective numbers at runtime.
 * NSGA2 with runtime objective numbers requires limited objNum.
 *
 * To change this value, define this macro before including this module, see the following code:
 *
 * \code {.cpp}
 * #define HEU_MAX_RUNTIME_OBJNUM 64
 * #include <unsupported/Eigen/CXX14/MetaHeuristic>
 * \endcode
 *
 */
#define HEU_MAX_RUNTIME_OBJNUM 32
#endif

namespace heu {
namespace internal {

/**
 * \ingroup HEU_GENETIC
 * \class MOGABase
 * \brief Base class for multiple objective genetic solvers.
 *
 * MOGABase maintains the numbers of objectives.
 *
 * This step of inheritance aims to support solvers with fixed and dynamic objective numbers. This template class has a
 * default implementation for fixed objective numbers and a partial specialization for dynamic. Lots of code-copying can
 * be avoided by such inheritance and parital specialization.
 *
 * \tparam Var_t
 * \tparam ObjNum
 * \tparam fOpt
 * \tparam rOpt
 * \tparam Args_t
 * \tparam _iFun_
 * \tparam _fFun_
 * \tparam _cFun_
 * \tparam _mFun_
 */
template <typename Var_t, int ObjNum, FitnessOption fOpt, RecordOption rOpt, class Args_t,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::initializeFun _iFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::fitnessFun _fFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::crossoverFun _cFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::mutateFun _mFun_>
class MOGABase : public MOGAAbstract<Var_t, ObjNum, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
  using Base_t = MOGAAbstract<Var_t, ObjNum, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

 public:
  ~MOGABase() {}
  HEU_MAKE_GABASE_TYPES(Base_t)

  /**
   * \brief Returns numbers of objectives
   *
   * \return constexpr size_t Numbers of objecvites at compile time.
   */
  constexpr size_t objectiveNum() const { return ObjNum; }
};

/**
 * \ingroup HEU_GENETIC
 * \class MOGABase
 * \brief Base class for multiple objective genetic solvers.
 *
 * \tparam Var_t
 * \tparam fOpt
 * \tparam rOpt
 * \tparam Args_t
 * \tparam _iFun_
 * \tparam _fFun_
 * \tparam _cFun_
 * \tparam _mFun_
 */
template <typename Var_t, FitnessOption fOpt, RecordOption rOpt, class Args_t,
          typename GAAbstract<Var_t, Eigen::ArrayXd, Args_t>::initializeFun _iFun_,
          typename GAAbstract<Var_t, Eigen::ArrayXd, Args_t>::fitnessFun _fFun_,
          typename GAAbstract<Var_t, Eigen::ArrayXd, Args_t>::crossoverFun _cFun_,
          typename GAAbstract<Var_t, Eigen::ArrayXd, Args_t>::mutateFun _mFun_>
class MOGABase<Var_t, Eigen::Dynamic, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>
    : public MOGAAbstract<Var_t, Eigen::Dynamic, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
 public:
  using Base_t = MOGAAbstract<Var_t, Eigen::Dynamic, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;
  HEU_MAKE_GABASE_TYPES(Base_t)

  MOGABase() { _objectiveNum = 0; }
  ~MOGABase() {}

  /**
   * \brief Return the value of objecvites
   *
   * \return int Number of objectives.
   */
  inline int objectiveNum() const { return _objectiveNum; }

  /**
   * \brief Set the Objective Num object
   *
   * \note Runtime assertion will fail if given _objNum is less than 2, or it exceeds HEU_MAX_RUNTIME_OBJNUM.
   * \note This member function exists only when template parameter ObjNum is Eigen::Dynamic.
   *
   * \param _objNum Number of objectives
   */
  inline void setObjectiveNum(int _objNum) {
    assert(_objNum > 1);
    assert(_objNum <= HEU_MAX_RUNTIME_OBJNUM);
    _objectiveNum = _objNum;
  }

 protected:
  int _objectiveNum;  ///< Default value is 0
};

}  //  namespace internal

}  //  namespace heu

#endif  // HEU_ MOGABASE_HPP
