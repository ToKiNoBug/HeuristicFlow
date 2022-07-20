/*
 Copyright © 2021-2022  TokiNoBug
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

#ifndef HEU_PSOABSTRCAT_HPP
#define HEU_PSOABSTRCAT_HPP

#ifdef HEU_DO_OUTPUT
#include <iostream>
#endif

#include "InternalHeaderCheck.h"
#include <HeuristicFlow/Global>
#include "PSOOption.hpp"
#include "PSOParameterPack.hpp"
#include "BoxWithVelocity.hpp"

namespace heu {

namespace internal {

/**
 * \ingroup HEU_PSO
 * \class PSOAbstract
 * \brief Internal base class for PSO solvers.
 *
 * Some fundamental types and functions are here.
 *
 * This template class has a specilization for PSO solvers without recording.
 *
 * \tparam Var_t Type of determination vector.
 * \tparam Fitness_t Type of fitness value.
 * \tparam Record Record trainning curve or not.
 * \tparam Arg_t Any other parameters
 * \tparam _iFun_ Initialization function at compile time
 * \tparam _fFun_ Initialization function at compile time
 *
 * \sa GABase It's counterpart in Genetic module.
 */
template <class Var_t, class Fitness_t, RecordOption Record, class Arg_t, BoxShape BS,
          typename PSOParameterPack<Var_t, Fitness_t, Arg_t>::iFun_t _iFun_,
          typename PSOParameterPack<Var_t, Fitness_t, Arg_t>::fFun_t _fFun_>
class PSOAbstract : public PSOParameterPack<Var_t, Fitness_t, Arg_t>,
                    public PSOParameterPack<Var_t, Fitness_t, Arg_t>::template iFunBody<_iFun_>,
                    public PSOParameterPack<Var_t, Fitness_t, Arg_t>::template fFunBody<_fFun_>,
                    public Box4PSO<Var_t, BS> {
  using Base_t = PSOParameterPack<Var_t, Fitness_t, Arg_t>;

 public:
  ~PSOAbstract() = default;
  HEU_MAKE_PSOPARAMETERPACK_TYPES(Base_t)
  using Scalar_t = typename toElement<Var_t>::type;

  /**
   * \brief Point is a pair of position together with fitness. It's speedless and.
   *
   */
  struct Point {
   public:
    /// The position of a point
    Var_t position;
    /// The fitness value of a point
    Fitness_t fitness;
  };

  /**
   * \brief Particle is a moveable point, it has its velocity and knows the best point it has ever
   * reached.
   *
   * \note Particle is inherited from Point, since it's just a point with speed and pBst.
   */
  struct Particle : public Point {
   public:
    /// The velocity of a point
    Var_t velocity;
    /// The ever reached best point.
    Point pBest;
  };

 public:
  /**
   * \brief Set the option object
   *
   * \param opt Option of PSO solver
   */
  inline void setOption(const PSOOption& opt) noexcept { _option = opt; }

  /**
   * \brief Get the option object
   *
   * \return const PSOOption& A const-ref to the option object
   */
  inline const PSOOption& option() const noexcept { return _option; }

  /**
   * \brief Get the generation.
   *
   * \return size_t generation
   */
  inline size_t generation() const noexcept { return _generation; }

  /**
   * \brief Get the fail times.
   *
   * Fail times refers to continuous generations that the solver failed to find a better solution.
   *
   * \return size_t fail times.
   */
  inline size_t failTimes() const noexcept { return _failTimes; }

  /**
   * \brief Get the population
   *
   * \return const std::vector<Particle>& A constant reference to the population
   */
  inline const std::vector<Particle>& population() const noexcept { return _population; }

  /**
   * \brief Get the global best solution that PSO has ever found.
   *
   * \return const Point& A const-ref to gBest.
   */
  inline const Point& globalBest() const noexcept { return gBest; }

  /*
   * \brief Set the range of position and velocity
   *
   * \param pMin Min val of position
   * \param pMax Max val of position
   * \param vMax Max val of velocity

  inline void setPVRange(const Var_t& pMin, const Var_t& pMax, const Var_t& vMax) noexcept {
    _posMin = pMin;
    _posMax = pMax;
    _velocityMax = vMax;
  }*/

  /*
   * \brief Set the range of position and velocity.
   *
   * \note This function will shape the box to a square box. Non't call this if you need a
   * non-square box.
   *
   * \param pMin Minium position value
   * \param pMax Maximum position value
   * \param vMax Maximum velocity absolute value

  inline void setPVRange(Scalar_t pMin, Scalar_t pMax, Scalar_t vMax) noexcept {
    for (int i = 0; i < this->_posMin.size(); i++) {
      this->_posMin[i] = pMin;
      this->_posMax[i] = pMax;
      this->_velocityMax[i] = vMax;
    }
  }*/

  /**
   * \brief Initialize the whole population
   *
   * This function reset the generation and failTimes to 0, and gBest to the first member of
   * population.
   *
   */
  void initializePop() noexcept {
    _population.resize(_option.populationSize);

    for (Particle& i : _population) {
      PSOExecutor<Base_t::HasParameters>::doInitialize(this, &i.position, &i.velocity);

      PSOExecutor<Base_t::HasParameters>::doFitness(this, &i.position, &i.fitness);

      i.pBest = i;
    }

    gBest = _population.front();
    _generation = 0;
    _failTimes = 0;
  }

 protected:
  /// The option of PSO solver
  PSOOption _option;

  /// Generation used.
  size_t _generation;

  /// failtimes
  size_t _failTimes;

  /*
/// Minimum position
Var_t _posMin;

/// Maximum position
Var_t _posMax;

/// Maximum velocity (absolute value)
Var_t _velocityMax;
*/

  /// All partiles in a vector
  std::vector<Particle> _population;

  /// The global pBest that the solver has ever found
  Point gBest;

  /**
   * \brief run the algorithm
   *
   * \tparam this_t Type of solver. This can be a PSOAbstract, or the solver type at the end of the
   * inheriting chain.
   *
   * __impl_run is designed to be a template function inorder to achieve compile polymorphism, kind
   * of like CRTP
   *
   * \sa GABase::run
   */
  template <class this_t = PSOAbstract>
  void __impl_run() noexcept {
    _generation = 0;
    _failTimes = 0;

    static_cast<this_t*>(this)->__impl_clearRecord();

    while (true) {
      _generation++;
      static_cast<this_t*>(this)->__impl_computeAllFitness();
      static_cast<this_t*>(this)->__impl_updatePGBest();

      static_cast<this_t*>(this)->template __impl_recordFitness<this_t>();
      if (_generation > _option.maxGeneration) {
#ifdef HEU_DO_OUTPUT
        std::cout << "Terminated by max generation limit" << std::endl;
#endif
        break;
      }

      if (_option.maxFailTimes > 0 && _failTimes > _option.maxFailTimes) {
#ifdef HEU_DO_OUTPUT
        std::cout << "Terminated by max failTime limit" << std::endl;
#endif
        break;
      }
#ifdef HEU_DO_OUTPUT
      std::cout << "Generation "
                << _generation
                //<<" , elite fitness="<<_eliteIt->fitness()
                << std::endl;
#endif
      static_cast<this_t*>(this)->__impl_updatePopulation();
    }
    _generation--;
  }

  /**
   * \brief Record fitness for non-recording solvers.
   * This function is useless here but it will be reloaded for PSOAbstract with recording.
   */
  inline void __impl_clearRecord() noexcept {}

  /**
   * \brief Record fitness for non-recording solvers.
   * This function is useless here but it will be reloaded for PSOAbstract with recording.
   */
  template <class this_t>
  inline void __impl_recordFitness() noexcept {}

  /**
   * \brief Compute fitness for the whole population
   *
   * In default cases, this function will boost the fitness computation via multi-threading.
   */
  void __impl_computeAllFitness() noexcept {
#ifdef HEU_HAS_OPENMP
    static const int32_t thN = threadNum();
#pragma omp parallel for schedule(dynamic, _population.size() / thN)
    for (int i = 0; i < _population.size(); i++) {
      Particle* ptr = &_population[i];
      PSOExecutor<Base_t::HasParameters>::doFitness(this, &ptr->position, &ptr->fitness);
    }
#else
    for (Particle& i : _population) {
      PSOExecutor<Base_t::HasParameters>::doFitness(this, &i.position, &i.fitness);
    }
#endif
  }

  // reloaded by template parameters to fit all types of `Args_t`
  template <bool _HasParameters, class unused = void>
  struct PSOExecutor {
    inline static void doInitialize(PSOAbstract* s, Var_t* pos, Var_t* velocity) noexcept {
      if constexpr (_iFun_ ==
                    PSOParameterPack<Var_t, Fitness_t,
                                     Arg_t>::defaultInitializeFunctionThatShouldNotBeCalled) {
        PSOExecutor<false, unused>::doInitiailize(s, pos, velocity);
      } else {
        s->runiFun(pos, velocity, &s->_arg);
      }
    }

    inline static void doFitness(PSOAbstract* s, const Var_t* pos, Fitness_t* f) noexcept {
      s->runfFun(pos, &s->_arg, f);
    }

    static_assert(Base_t::HasParameters == _HasParameters,
                  "A wrong specialization of PSOExecuter is called");
  };

  template <class unused>
  struct PSOExecutor<false, unused> {
    inline static void doInitialize(PSOAbstract* s, Var_t* pos, Var_t* velocity) noexcept {
      if constexpr (_iFun_ ==
                    PSOParameterPack<Var_t, Fitness_t,
                                     Arg_t>::defaultInitializeFunctionThatShouldNotBeCalled) {
        if constexpr (!array_traits<Var_t>::isFixedSize) {
          if constexpr (array_traits<Var_t>::isEigenClass) {
            // resize for Eigen matrices
            if constexpr (array_traits<Var_t>::isMatrix) {
              pos->setZero(s->boxRows(), s->boxCols());
              velocity->setZero(s->boxRows(), s->boxCols());
            } else if constexpr (array_traits<Var_t>::isColVector) {
              pos->setZero(s->dimensions(), 1);
              velocity->setZero(s->dimensions(), 1);
            } else {  // must be row vector
              pos->setZero(1, s->dimensions());
              velocity->setZero(1, s->dimensions());
            }

          } else {
            // resize for std vectors
            pos->resize(s->dimensions());
            velocity->resize(s->dimensions());
          }
        }
        // evulate
        for (int idx = 0; idx < s->dimensions(); idx++) {
          at(*velocity, idx) = 0;
          if constexpr (BS == BoxShape::SQUARE_BOX) {
            at(*pos, idx) = randD(s->posMin(), s->posMax());

          } else {
            at(*pos, idx) = randD(s->posMin()[idx], s->posMax()[idx]);
          }
        }
      } else {
        s->runiFun(pos, velocity);
      }
    }

    inline static void doFitness(PSOAbstract* s, const Var_t* pos, Fitness_t* f) noexcept {
      s->runfFun(pos, f);
    }

    static_assert(Base_t::HasParameters == false,
                  "A wrong specialization of PSOExecuter is called");
  };
};

#define HEU_MAKE_PSOABSTRACT_TYPES(Base_t)    \
  HEU_MAKE_PSOPARAMETERPACK_TYPES(Base_t)     \
  using Point_t = typename Base_t::Point;     \
  using Scalar_t = typename Base_t::Scalar_t; \
  using Particle_t = typename Base_t::Particle;

/**
 * \ingroup HEU_PSO
 * \class PSOAbstract<Var_t, Fitness_t, RECORD_FITNESS, Arg_t, _iFun_, _fFun_>
 * \brief partial specialization for PSO with recording
 *
 * \note PSOAbstract with record is herited from PSOAbstract without record.
 *
 * \tparam Var_t Type of determination vector.
 * \tparam Fitness_t Type of fitness value.
 * \tparam Arg_t Any other parameters
 * \tparam _iFun_ Initialization function at compile time
 * \tparam _fFun_ Initialization function at compile time
 */
template <class Var_t, class Fitness_t, class Arg_t, BoxShape BS,
          typename PSOParameterPack<Var_t, Fitness_t, Arg_t>::iFun_t _iFun_,
          typename PSOParameterPack<Var_t, Fitness_t, Arg_t>::fFun_t _fFun_>
class PSOAbstract<Var_t, Fitness_t, RECORD_FITNESS, Arg_t, BS, _iFun_, _fFun_>
    : public PSOAbstract<Var_t, Fitness_t, DONT_RECORD_FITNESS, Arg_t, BS, _iFun_, _fFun_> {
  using Base_t = PSOAbstract<Var_t, Fitness_t, DONT_RECORD_FITNESS, Arg_t, BS, _iFun_, _fFun_>;
  friend Base_t;

 public:
  ~PSOAbstract() = default;
  HEU_MAKE_PSOABSTRACT_TYPES(Base_t)

  /**
   * \brief Get the fitness record
   *
   * \return const std::vector<Fitness_t>& The fitness record.
   */
  const std::vector<Fitness_t>& record() const noexcept { return _record; }

 protected:
  /// The fitness record
  std::vector<Fitness_t> _record;

  /**
   * \brief Clear the record.
   *
   * \sa PSOAbstract::__impl_clearRecord
   *
   */
  inline void __impl_clearRecord() noexcept {
    _record.clear();
    _record.reserve(this->_option.maxGeneration + 1);
  }

  /**
   * \brief Record fitness.
   *
   * \sa PSOAbstract::__impl_recordFitness
   *
   */
  template <class this_t>
  inline void __impl_recordFitness() noexcept {
    _record.emplace_back(static_cast<this_t*>(this)->bestFitness());
  }
};

}  //  namespace internal

}  //  namespace heu

#endif  // HEU_PSOABSTRCAT_HPP
