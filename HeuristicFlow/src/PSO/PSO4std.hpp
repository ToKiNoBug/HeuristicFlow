#ifndef HEU_PSO4STD_HPP
#define HEU_PSO4STD_HPP

#include <array>
#include <vector>
#include <tuple>
#include <type_traits>

#include "InternalHeaderCheck.h"
#include "PSOOption.hpp"
#include "PSOBase.hpp"

namespace heu {

namespace internal {

template <typename Var_t, int DIM, FitnessOption FitnessOpt, RecordOption RecordOpt, class Arg_t,
          typename internal::PSOParameterPack<Var_t, double, Arg_t>::iFun_t _iFun_,
          typename internal::PSOParameterPack<Var_t, double, Arg_t>::fFun_t _fFun_>
class PSO4std : public internal::PSOBase<Var_t, DIM, double, RecordOpt, Arg_t, _iFun_, _fFun_> {
  using Base_t = internal::PSOBase<Var_t, DIM, double, RecordOpt, Arg_t, _iFun_, _fFun_>;

 public:
  PSO4std() {}
  ~PSO4std() {}
  HEU_MAKE_PSOABSTRACT_TYPES(Base_t)
  friend class internal::PSOAbstract<Var_t, double, DONT_RECORD_FITNESS, Arg_t, _iFun_, _fFun_>;

  static constexpr ContainerOption Flag =
      (std::is_same<Var_t, stdVecD_t<DIM>>::value) ? ContainerOption::Std : ContainerOption::Custom;

 protected:
  /**
   * \brief To determine whether fitness a is better than b
   *
   * \param a The first input
   * \param b The second input
   * \return true A is better than B
   * \return false A is not better than B
   */
  inline static bool isBetterThan(double a, double b) {
    if (FitnessOpt == FitnessOption::FITNESS_GREATER_BETTER) {
      return a > b;
    } else {
      return a < b;
    }
  }

  /**
   * \brief Update gBest and pBest
   *
   */
  void __impl_updatePGBest() {
    // gBest for current generation
    Point_t* curGBest = &this->_population.front().pBest;

    for (Particle_t& i : this->_population) {
      if (isBetterThan(i.fitness, i.pBest.fitness)) {
        i.pBest = i;
      }

      if (isBetterThan(i.pBest.fitness, curGBest->fitness)) {
        curGBest = &i.pBest;
      }
    }

    if (isBetterThan(curGBest->fitness, this->gBest.fitness)) {
      this->_failTimes = 0;
      this->gBest = *curGBest;
    } else {
      this->_failTimes++;
    }
  }

  /**
   * \brief Update the position and velocity of all particles
   *
   */
  void __impl_updatePopulation() {
#ifdef EIGEN_HAS_OPENMP
    static const int32_t thN = Eigen::nbThreads();
#pragma omp parallel for schedule(dynamic, this->_population.size() / thN)
#endif  //  EIGEN_HAS_OPENMP
    for (int index = 0; index < this->_population.size(); index++) {
      Particle_t& i = this->_population[index];
      const double rndP = randD();
      const double rndG = randD();
      for (int idx = 0; idx < this->dimensions(); idx++) {
        i.velocity[idx] = this->_option.inertiaFactor * i.velocity[idx] +
                          this->_option.learnFactorP * rndP * (i.pBest.position[idx] - i.position[idx]) +
                          this->_option.learnFactorG * rndG * (this->gBest.position[idx] - i.position[idx]);
        if (std::abs(i.velocity[idx]) > this->_velocityMax[idx]) {
          i.velocity[idx] = sign(i.velocity[idx]) * this->_velocityMax[idx];
        }
        i.position[idx] += i.velocity[idx];
        i.position[idx] = std::max(i.position[idx], this->_posMin[idx]);
        i.position[idx] = std::min(i.position[idx], this->_posMax[idx]);
      }
    }
  }

 private:
  static_assert(!(std::is_scalar<Var_t>::value), "Var_t should be a non-scalar type");
  static_assert(DIM != 0, "Template parameter DIM cannot be 0. For dynamic dims, use Eigen::Dynamic");
  static_assert(DIM > 0 || DIM == Eigen::Dynamic, "Invalid template parameter DIM");
  // static_assert(isEigenTypes == false, "Wrong specialization of PSO");
};

}  //  namespace internal

}  //  namespace heu

#endif  //  HEU_PSO4STD_HPP
