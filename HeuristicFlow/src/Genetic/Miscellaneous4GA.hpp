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

#ifndef HEU_MISCELLANEOUS4GA_HPP
#define HEU_MISCELLANEOUS4GA_HPP

#include <assert.h>
#include <HeuristicFlow/Global>
#include "InternalHeaderCheck.h"
#include "GAAbstract.hpp"

namespace heu {

namespace internal {

template <typename Var_t, ContainerOption dvo>
struct imp_GADefaults_DVO {
  template <DivCode _r>
  inline static void imp_cFunNd(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2) noexcept {
    static const double constexpr r = DivDecode<_r>::real;
    static_assert(r > 0, "r shouldn't be less than 0");
    static_assert(r < 1, "r shouldn't be greater than 1");

    const size_t N = p1->size();
    for (size_t i = 0; i < N; i++) {
      c1->operator[](i) = r * (p1->operator[](i)) + (1 - r) * (p2->operator[](i));
      c2->operator[](i) = r * (p2->operator[](i)) + (1 - r) * (p1->operator[](i));
    }
  }

  inline static void imp_cFunSwapNs(const Var_t *p1, const Var_t *p2, Var_t *c1,
                                    Var_t *c2) noexcept {
    const size_t N = p1->size();
    const size_t idx = randIdx(N);

    c1->topRows(idx) = p1->topRows(idx);
    c2->topRows(idx) = p2->topRows(idx);
    c1->bottomRows(N - idx) = p2->bottomRows(N - idx);
    c2->bottomRows(N - idx) = p1->bottomRows(N - idx);
  }
};

template <typename Var_t>
struct imp_GADefaults_DVO<Var_t, ContainerOption::Eigen> {
  template <DivCode _r>
  inline static void imp_cFunNd(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2) noexcept {
    static const double constexpr r = DivDecode<_r>::real;
    static_assert(r > 0, "r shouldn't be less than 0");
    static_assert(r < 1, "r shouldn't be greater than 1");

    *c1 = r * (*p1) + (1 - r) * (*p2);
    *c2 = r * (*p2) + (1 - r) * (*p1);
  }

  inline static void imp_cFunSwapNs(const Var_t *p1, const Var_t *p2, Var_t *c1,
                                    Var_t *c2) noexcept {
    const size_t N = p1->size();
    const size_t idx = randIdx(N);

    for (size_t i = 0; i < N; i++) {
      if (i < idx) {
        c1->operator[](i) = p1->operator[](i);
        c2->operator[](i) = p2->operator[](i);
      } else {
        c1->operator[](i) = p2->operator[](i);
        c2->operator[](i) = p1->operator[](i);
      }
    }
  }
};

}  // namespace internal

/**
 * \ingroup HEU_GENETIC
 * \struct GADefaults
 * \brief The GADefaults struct defines several candidate operations for GA.
 *
 * \tparam Var_t type of determinate vector
 * \tparam Args_t type of other parameters in genetic solver (void as default)
 * \tparam dvo double vector option of Var_t (cpp std as default)
 *
 * The format of function starts with (i/c/m)Fun(N/X/_)(d/b/s).
 * - In the first brace, i means initialization, while c for crossver and m for mutation.
 * - In the second brace, N means fix-sized vectors and X for dynamic sized vectors. While _ means
 * both are supported.
 * - In the last brace, d means `double`(real encoding), wile b for `bool` (binaray encoding) and s
 * for symbolic encoding.
 *
 * There isn't function for fFun because that it right the problem you will solve.
 *
 * \note `cFunNd` and `cFunXd` are only avaliable for real encoding, while the rest crossover
 * function supports all encodings.
 *
 * \note All iFun and mFun requires `Args_t` to be (or inherit from) a box-constraint type.
 * Otherwise static assertion will fail.
 *
 * \note To avoid static assertion failuer with non-box constraint types of `Args_t`, some functions
 * have a unused parameter to stop unnecessary template instatiation. For example, when using
 * `iFunNd`, use `iFunNd<>`. Don't miss the template braces if this static member function is a
 * templated one.
 *
 * \note This struct has a specilization for `Args_t` as `void` (which means no parameters).
 *
 * \sa GADefaults<Var_t, void, dvo>
 * \sa internal::GABase internal::GAAbstract
 *
 */
template <typename Var_t, class Args_t = void, ContainerOption dvo = ContainerOption::Std>
struct GADefaults {
 private:
  static_assert(!std::is_same<Args_t, void>::value,
                "The compiler run into a incorrect branch of partial specialization");

  template <BoxShape BS, typename unused = void>
  struct RealBoxOp {
    static_assert(Args_t::Shape == BoxShape::RECTANGLE_BOX, "Wrong specialization");
    // non-square box

   private:
  };

  template <typename unused>
  struct RealBoxOp<BoxShape::SQUARE_BOX, unused> {
    // square box
  };

  template <BoxShape BS, typename unused = void>
  struct SymBoxOp  //  non-square box
  {
    using scalar_t = typename toElement<Var_t>::type;

    inline static void imp_domFunNs(const Var_t *src, Var_t *v, const Args_t *box) noexcept {
      *v = *src;
      const size_t idx = randIdx(v->size());
      const auto val = v->operator[](idx);
      const size_t numLess = val - box->min()[idx];
      const size_t numGreater = box->max()[idx] - val;

      if (randD() * (numLess + numGreater) <= numLess)
        v->operator[](idx) = randIdx(box->min()[idx], val);
      else
        v->operator[](idx) = randIdx(scalar_t(val + 1), scalar_t(box->max()[idx] + 1));
    }
  };

  template <typename unused>
  struct SymBoxOp<BoxShape::SQUARE_BOX, unused> {
    using scalar_t = typename toElement<Var_t>::type;

    inline static void imp_domFunNs(const Var_t *src, Var_t *v, const Args_t *box) noexcept {
      *v = *src;
      const size_t idx = randIdx(v->size());
      const auto val = v->operator[](idx);
      const size_t numLess = val - box->min();
      const size_t numGreater = box->max() - val;

      if (randD() * (numLess + numGreater) <= numLess)
        v->operator[](idx) = randIdx(box->min(), val);
      else
        v->operator[](idx) = randIdx(scalar_t(val + 1), scalar_t(box->max() + 1));
    }
  };

 public:
  /**
   * \brief Univerial initialize function.
   *
   * \note This function requires `Args_t` to be (or derived from) a box-constraint.
   *
   * \param v pointer of Var_t
   * \param box Pointer to the box constraint
   *
   */
  inline static void iFun(Var_t *v, const Args_t *box) noexcept { box->initialize(v); }

  /**
   * \brief Default crossover function for fixed-size float/double array/vector
   *
   * Real encoding enables arithmatical crossver like this:\n
   * c1=r*p1+(1-r)*p2 and c2=r*p2+(1-r)*p1.\n
   * Where r is a ratio in range [0,1]. The less r is, the more similarity c1 and p1 are.
   *
   * \tparam _r Crossover ratio r, encoded as DivCode. Default value is 0.2.
   * \param p1 The first parent
   * \param p2 The second parnet
   * \param c1 The first child
   * \param c2 The second child
   */
  template <DivCode _r = DivEncode<1, 5>::code>
  inline static void cFunNd(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2,
                            const Args_t *) noexcept {
    GADefaults<Var_t, void, dvo>::template cFunNd<_r>(p1, p2, c1, c2);
  }

  /**
   * \brief Default crossover function for dynamic-size float/double array/vector
   *
   * Real encoding enables arithmatical crossver like this:\n
   * c1=r*p1+(1-r)*p2 and c2=r*p2+(1-r)*p1.\n
   * Where r is a ratio in range [0,1]. The less r is, the more similarity c1 and p1 are.
   *
   * \tparam _r Crossover ratio r, encoded as DivCode. Default value is 0.2.
   * \param p1 The first parent
   * \param p2 The second parnet
   * \param c1 The first child
   * \param c2 The second child
   */
  template <DivCode _r = DivEncode<1, 5>::code>
  inline static void cFunXd(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2,
                            const Args_t *) noexcept {
    GADefaults<Var_t, void, dvo>::template cFunXd<_r>(p1, p2, c1, c2);
  }

  /**
   * \brief Default crossover function for fixed-size array/vector (Genetic with args)
   *
   * This is the original crossover method suitable for binary, symbolic and real vectors. It acts
   * like chromosomes.
   *
   * \param p1 The first parent
   * \param p2 The second parnet
   * \param c1 The first child
   * \param c2 The second child
   */
  inline static void cFunSwapNs(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2,
                                const Args_t *) noexcept {
    GADefaults<Var_t, void, dvo>::cFunSwapNs(p1, p2, c1, c2);
  }

  /**
   * \brief Default crossover function for runtime-size array/vector (Genetic with args)
   *
   * This is the original crossover method suitable for binary, symbolic and real vectors. It acts
   * like chromosomes.
   *
   * \param p1 The first parent
   * \param p2 The second parnet
   * \param c1 The first child
   * \param c2 The second child
   */
  inline static void cFunSwapXs(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2,
                                const Args_t *) noexcept {
    GADefaults<Var_t, void, dvo>::cFunSwapXs(p1, p2, c1, c2);
  }

  /**
   * \brief Discrete random selection crossover by probability for fixed-size array/vector (with
   * args)
   *
   * Child's value of each element is chosen from their parents stochastically by probability p.
   * Where p is the probability that c1 choose its value from p1 and c2 from p2.
   *
   * This crossover method is suitable for all encoding types.
   *
   * \tparam p Encoded in DivCode. Default value 0.5
   * \param p1 The first parent
   * \param p2 The second parnet
   * \param c1 The first child
   * \param c2 The second child
   */
  template <DivCode p = DivCode::DivCode_Half>
  inline static void cFunRandNs(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2,
                                const Args_t *) noexcept {
    GADefaults<Var_t, void, dvo>::template cFunRandNs<p>(p1, p2, c1, c2);
  }

  /**
   * \brief Discrete random selection crossover by probability for fixed-size array/vector (with
   * args)
   *
   * Child's value of each element is chosen from their parents stochastically by probability p.
   * Where p is the probability that c1 choose its value from p1 and c2 from p2.
   *
   * This crossover method is suitable for all encoding types.
   *
   * \tparam posCode Encoded in DivCode. Default value 0.5
   * \param p1 The first parent
   * \param p2 The second parnet
   * \param c1 The first child
   * \param c2 The second child
   */
  template <DivCode posCode = DivCode::DivCode_Half>
  inline static void cFunRandXs(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2,
                                const Args_t *) noexcept {
    GADefaults<Var_t, void, dvo>::template cFunRandXs<posCode>(p1, p2, c1, c2);
  }

  inline static void mFun(const Var_t *src, Var_t *v, const Args_t *box) noexcept {
    *v = *src;
    box->applyDelta(v);
  }

  /**
   * \brief Default mutate function for symbolic vectors (fixed and runtime size)
   *
   * \tparam unused unused template parameter is introduced to avoid static assertion failuer when
   * non-box Args_t is used.
   *
   * \note This function requires `Args_t` to be (or derived from) a symbolic box-constraint.
   *
   * \param src The gene before mutation
   * \param v The gene after mutation
   * \param box The box constraint
   */
  template <typename unused = void>
  inline static void mFun_s(const Var_t *src, Var_t *v, const Args_t *box) noexcept {
    static_assert(Args_t::isBox, "Default mFun requires Args_t to be a box constriant");
    static_assert(Args_t::Encoding == EncodeType::Symbolic, "mFun_s requires symbolic box");
    static_assert(std::is_same<typename Args_t::Var_t, Var_t>::value,
                  "Box and Var_t types must be same");

    SymBoxOp<Args_t::Shape>::imp_domFunNs(src, v, box);
  }
};

template <typename Var_t, ContainerOption dvo>
struct GADefaults<Var_t, void, dvo> {
  /**
   * \brief Initialization function for fixed double array without args.
   *
   * iFun without args as box constraint can't provide much support for initialization.
   * In this function, the minimum and maximum can only be passed through template parameters.
   *
   * This function works like initialization with a square real box.
   *
   * \tparam _min Minimum value in DivCode (0.0 for default)
   * \tparam _max Maximum value in DivCode (1.0 for default)
   * \param p Var_t to be initialized.
   */
  template <DivCode _min = DivEncode<0, 1>::code, DivCode _max = DivEncode<1, 1>::code>
  inline static void iFunNd(Var_t *p) noexcept {
    static constexpr double min = DivDecode<_min>::real;
    static constexpr double max = DivDecode<_max>::real;
    // static const constexpr bool isValid=(max>min);
    // static_assert(isValid,"Max should be greater than min");

    for (int64_t idx = 0; idx < p->size(); idx++) {
      p->operator[](idx) = randD(min, max);
    }
  }

  /**
   * \brief Initialization function for fixed float array without args.
   *
   * iFun without args as box constraint can't provide much support for initialization.
   * In this function, the minimum and maximum can only be passed through template parameters.
   *
   * This function works like initialization with a square real box.
   *
   * \tparam _min Minimum value in DivCode (0.0 for default)
   * \tparam _max Maximum value in DivCode (1.0 for default)
   * \param p Var_t to be initialized.
   */
  template <DivCode _min = DivEncode<0, 1>::code, DivCode _max = DivEncode<1, 1>::code>
  inline static void iFunNf(Var_t *p) noexcept {
    iFunNd<_min, _max>(p);
  }

  /**
   * \brief This function has the same effect as its counterpart in GADefaults.
   *
   * \sa GADefaults::cFunNd
   */
  template <DivCode _r = DivEncode<1, 5>::code>
  inline static void cFunNd(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2) noexcept {
    internal::template imp_GADefaults_DVO<Var_t, dvo>::template imp_cFunNd<_r>(p1, p2, c1, c2);
  }

  /**
   * \brief This function has the same effect as its counterpart in GADefaults.
   *
   * \sa GADefaults::cFunXd
   */
  template <DivCode _r = DivEncode<1, 5>::code>
  inline static void cFunXd(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2) noexcept {
#define Heu_PRIVATE_IMP_cFunX \
  c1->resize(p1->size());     \
  c2->resize(p2->size());

    Heu_PRIVATE_IMP_cFunX

        cFunNd<_r>(p1, p2, c1, c2);
  }

  /**
   * \brief This function has the same effect as its counterpart in GADefaults.
   *
   * \sa GADefaults::cFunSwapNs
   */
  inline static void cFunSwapNs(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2) noexcept {
    internal::template imp_GADefaults_DVO<Var_t, dvo>::imp_cFunSwapNs(p1, p2, c1, c2);
  }

  /**
   * \brief This function has the same effect as its counterpart in GADefaults.
   *
   * \sa GADefaults::cFunSwapXs
   */
  inline static void cFunSwapXs(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2) noexcept {
    Heu_PRIVATE_IMP_cFunX

        cFunSwapNs(p1, p2, c1, c2);
  }

  /**
   * \brief This function has the same effect as its counterpart in GADefaults.
   *
   * \sa GADefaults::cFunRandNs
   */
  template <DivCode p = DivCode::DivCode_Half>
  inline static void cFunRandNs(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2) noexcept {
    static const double constexpr r = DivDecode<p>::real;
    static_assert(r > 0, "A probability shoule be greater than 0");
    static_assert(r < 1, "A probability shoule be less than 1");
    const size_t N = p1->size();
    for (size_t i = 0; i < N; i++) {
      c1->operator[](i) = ((randD() < r) ? p1 : p2)->operator[](i);
      c2->operator[](i) = ((randD() < r) ? p2 : p1)->operator[](i);
    }
  }

  /**
   * \brief This function has the same effect as its counterpart in GADefaults.
   *
   * \sa GADefaults::cFunRandXs
   */
  template <DivCode p = DivCode::DivCode_Half>
  inline static void cFunRandXs(const Var_t *p1, const Var_t *p2, Var_t *c1, Var_t *c2) noexcept {
    Heu_PRIVATE_IMP_cFunX cFunRandNs(p1, p2, c1, c2);
  }
};

}  //  namespace heu

#endif  //  HEU_MISCELLANEOUS4GA_HPP
