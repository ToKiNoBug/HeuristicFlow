#ifndef HEU_BOXCONSTRAINT_HPP
#define HEU_BOXCONSTRAINT_HPP

#include <limits>
#include <type_traits>
#include "../../Global"

#include "SizeBody4BoxConstraint.hpp"

namespace heu {

namespace {
template <class Derived, class Var_t>
class BoxBase {
  // using Var_t = typename Derived::Var_t;
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;

 public:
  inline void initialize(Var_t* v) const noexcept {
    static_cast<const Derived*>(this)->initializeSize(v);
    for (int idx = 0; idx < v->size(); idx++) {
      if constexpr (std::is_same_v<Scalar_t, bool>) {
        at(*v, idx) = randIdx(2) == 0;
      } else if constexpr (std::is_integral_v<Scalar_t>) {
        at(*v, idx) = randIdx<Scalar_t>(static_cast<const Derived*>(this)->min(idx),
                                        1 + static_cast<const Derived*>(this)->max(idx));
      } else {
        at(*v, idx) = randD(static_cast<const Derived*>(this)->min(idx),
                            static_cast<const Derived*>(this)->max(idx));
      }
    }
  }

  inline void applyConstraint(Var_t* v) const noexcept {
    assert(v->size() == static_cast<const Derived*>(this)->dimensions());
    for (int idx = 0; idx < v->size(); idx++) {
      at(*v, idx) = std::min(at(*v, idx), static_cast<const Derived*>(this)->max(idx));
      at(*v, idx) = std::max(at(*v, idx), static_cast<const Derived*>(this)->min(idx));
    }
  }

  inline void applyDelta(Var_t* v) const noexcept {
    assert(v->size() == static_cast<const Derived*>(this)->dimensions());
    const int idx = randIdx(0, static_cast<const Derived*>(this)->dimensions());
    if constexpr (std::is_same_v<Scalar_t, bool>) {
      // binary mutation
      at(*v, idx) = !at(*v, idx);
    } else {
      // integer mutation
      at(*v, idx) += (randIdx(2) == 0) ? 1 : -1;
      at(*v, idx) = std::min(at(*v, idx), static_cast<const Derived*>(this)->max(idx));
      at(*v, idx) = std::max(at(*v, idx), static_cast<const Derived*>(this)->min(idx));
    }
  }

 protected:
  inline void assert4Size(const int idx) const noexcept {
    assert(idx >= 0 && idx < static_cast<const Derived*>(this)->dimensions());
  }
};

////////////////

template <class Derived, class Var_t>
class MatBoxBase {
 public:
 protected:
  inline void assert4Size(const int r, const int c) const noexcept {
    assert(r >= 0 && r < static_cast<const Derived*>(this)->boxRows());
    assert(c >= 0 && c < static_cast<const Derived*>(this)->boxCols());
  }
};

/////////////////////////////

template <class Var_t>
class SquareBoxCore {
 public:
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;

  inline Scalar_t& min() noexcept { return _minS; }
  inline Scalar_t min() const noexcept { return _minS; }
  inline Scalar_t min(const int) const noexcept { return _minS; }

  inline Scalar_t& max() noexcept { return _maxS; }
  inline Scalar_t max() const noexcept { return _maxS; }
  inline Scalar_t max(const int) const noexcept { return _maxS; }

  inline void setRange(const Scalar_t __min, const Scalar_t __max) noexcept {
    assert(__min < __max);
    _minS = __min;
    _maxS = __max;
  }

 private:
  Scalar_t _minS;
  Scalar_t _maxS;
};

template <class Var_t>
class BooleanBoxCore {
 public:
  static_assert(std::is_same_v<typename array_traits<Var_t>::Scalar_t, bool>);
  using Scalar_t = bool;

  inline constexpr bool min() const noexcept { return 0; }
  inline constexpr bool min(const int) const noexcept { return 0; }

  inline constexpr bool max() const noexcept { return 1; }
  inline constexpr bool max(const int) const noexcept { return 1; }
};

////////////////////////////

template <class Var>
class SquareBoxConstDimVec
    : public BoxBase<SquareBoxConstDimVec<Var>, Var>,
      public std::conditional_t<
          !std::is_same_v<bool, typename array_traits<Var>::Scalar_t>,
          std::conditional_t<!std::is_same_v<bool, typename array_traits<Var>::Scalar_t>,
                             SquareBoxCore<Var>, BooleanBoxCore<Var>>,
          BooleanBoxCore<Var>>,
      public internal::SquareBoxSizeBody<Var> {
 public:
  using Var_t = Var;
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;
  using BoxBase_t = BoxBase<SquareBoxConstDimVec<Var>, Var>;
  static_assert(array_traits<Var_t>::isFixedSize);
  static_assert(array_traits<Var_t>::isVector);
  /*
    inline constexpr int dimensions() const noexcept { return array_traits<Var_t>::sizeCT; }

    inline void initializeSize(Var_t*) const noexcept {}

   private:
   */
};

template <class Var>
class SquareBoxConstDimMat
    : public BoxBase<SquareBoxConstDimMat<Var>, Var>,
      public MatBoxBase<SquareBoxConstDimMat<Var>, Var>,
      public std::conditional_t<!std::is_same_v<bool, typename array_traits<Var>::Scalar_t>,
                                SquareBoxCore<Var>, BooleanBoxCore<Var>>,
      public internal::SquareBoxSizeBody<Var> {
 public:
  using Var_t = Var;
  using BoxBase_t = BoxBase<SquareBoxConstDimMat<Var>, Var>;
  using MatBoxBase_t = MatBoxBase<SquareBoxConstDimMat<Var>, Var>;
  static_assert(array_traits<Var_t>::isFixedSize);
  static_assert(!array_traits<Var_t>::isVector);

  using Scalar_t = typename array_traits<Var_t>::Scalar_t;
  /*
    inline constexpr int dimensions() const noexcept { return array_traits<Var_t>::sizeCT; }

    inline constexpr int boxRows() const noexcept { return array_traits<Var_t>::rowsCT; }

    inline constexpr int boxCols() const noexcept { return array_traits<Var_t>::colsCT; }

    inline void initializeSize(Var_t*) const noexcept {}
    */
};

template <class Var>
class SquareBoxDynamicDimVec
    : public BoxBase<SquareBoxDynamicDimVec<Var>, Var>,
      public std::conditional_t<!std::is_same_v<bool, typename array_traits<Var>::Scalar_t>,
                                SquareBoxCore<Var>, BooleanBoxCore<Var>>,
      public internal::SquareBoxSizeBody<Var> {
 public:
  using Var_t = Var;
  using BoxBase_t = BoxBase<SquareBoxConstDimVec<Var>, Var>;
  static_assert(!array_traits<Var_t>::isFixedSize);
  static_assert(array_traits<Var_t>::isVector);
  /*
    inline int dimensions() const noexcept { return _dimensions; }

    inline void setDimensions(const int dim) noexcept {
      assert(dim > 0);
      _dimensions = dim;
    }

    inline void initializeSize(Var_t* v) const noexcept {
      assert(_dimensions > 0);
      v->resize(_dimensions);
    }

   private:
    int _dimensions;*/
};

template <class Var>
class SquareBoxDynamicDimMat
    : public BoxBase<SquareBoxDynamicDimMat<Var>, Var>,
      public MatBoxBase<SquareBoxDynamicDimMat<Var>, Var>,
      public std::conditional_t<!std::is_same_v<bool, typename array_traits<Var>::Scalar_t>,
                                SquareBoxCore<Var>, BooleanBoxCore<Var>>,
      public internal::SquareBoxSizeBody<Var> {
 public:
  using Var_t = Var;
  using BoxBase_t = BoxBase<SquareBoxDynamicDimMat<Var>, Var>;
  using MatBoxBase_t = MatBoxBase<SquareBoxDynamicDimMat<Var>, Var>;
  static_assert(!array_traits<Var_t>::isFixedSize);
  static_assert(!array_traits<Var_t>::isVector);

  inline int dimensions() const noexcept { return _rows * _cols; }

  inline int boxRows() const noexcept { return _rows; }

  inline int boxCols() const noexcept { return _cols; }

  inline void setDimensions(const int rows, const int cols) noexcept {
    assert(rows > 0);
    assert(cols > 0);
    _rows = rows;
    _cols = cols;
  }

  inline void initializeSize(Var_t* v) const noexcept {
    assert(_rows > 0);
    assert(_cols > 0);
    v->resize(_rows, _cols);
  }

 private:
  int _rows;
  int _cols;
};

template <class Var_t>
using SquareBox_t = std::conditional_t<
    array_traits<Var_t>::isFixedSize,
    std::conditional_t<array_traits<Var_t>::isVector, SquareBoxConstDimVec<Var_t>,
                       SquareBoxConstDimMat<Var_t>>,

    std::conditional_t<array_traits<Var_t>::isVector, SquareBoxDynamicDimVec<Var_t>,
                       SquareBoxDynamicDimMat<Var_t>>>;

//////////////////////////////////////////

template <class Var_t>
class NonSquareBoxCore {
 public:
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;

  inline Var_t& min() noexcept { return _minV; }
  inline const Var_t& min() const noexcept { return _minV; }
  inline Scalar_t min(const int idx) const noexcept { return at(_minV, idx); }

  inline Var_t& max() noexcept { return _maxV; }
  inline const Var_t& max() const noexcept { return _maxV; }
  inline Scalar_t max(const int idx) const noexcept { return at(_maxV, idx); }

  inline void setRange(const Var_t& __min, const Var_t& __max) const noexcept {
    assert(__min.size() == __max.size());
    _minV = __min;
    _maxV = __max;
  }

 private:
  Var_t _minV;
  Var_t _maxV;
};
/////////////////

template <class Var>
class NonSquareBoxConstDimVec : public BoxBase<NonSquareBoxConstDimVec<Var>, Var>,
                                public NonSquareBoxCore<Var> {
 public:
  using Var_t = Var;
  using BoxBase_t = BoxBase<NonSquareBoxConstDimVec<Var>, Var>;
  static_assert(array_traits<Var_t>::isFixedSize);
  static_assert(array_traits<Var_t>::isVector);

  inline constexpr int dimensions() const noexcept { return array_traits<Var_t>::sizeCT; }

  inline void initializeSize(Var_t*) const noexcept {}
};

template <class Var>
class NonSquareBoxConstDimMat : public BoxBase<NonSquareBoxConstDimMat<Var>, Var>,
                                public MatBoxBase<NonSquareBoxConstDimMat<Var>, Var>,
                                public NonSquareBoxCore<Var> {
 public:
  using Var_t = Var;
  using BoxBase_t = BoxBase<NonSquareBoxConstDimMat<Var>, Var>;
  using MatBoxBase_t = MatBoxBase<NonSquareBoxConstDimMat<Var>, Var>;

  static_assert(array_traits<Var_t>::isFixedSize);
  static_assert(!array_traits<Var_t>::isVector);

  inline constexpr int dimensions() const noexcept { return array_traits<Var_t>::sizeCT; }

  inline constexpr int boxRows() const noexcept { return array_traits<Var_t>::rowsCT; }

  inline constexpr int boxCols() const noexcept { return array_traits<Var_t>::colsCT; }

  inline void initializeSize(Var_t* v) const noexcept {}
};

template <class Var>
class NonSquareBoxDynamicDimVec : public BoxBase<NonSquareBoxDynamicDimVec<Var>, Var>,
                                  public NonSquareBoxCore<Var> {
 public:
  using Var_t = Var;
  using BoxBase_t = BoxBase<NonSquareBoxDynamicDimVec<Var>, Var>;
  static_assert(!array_traits<Var_t>::isFixedSize);
  static_assert(array_traits<Var_t>::isVector);

  inline int dimensions() const noexcept {
    assert(this->min().size() == this->max().size());
    return this->min().size();
  }

  inline int setDimensions(const int dim) const noexcept {
    assert(dim > 0);
    this->min().resize(dim);
    this->max().resize(dim);
  }

  inline void initializeSize(Var_t* v) const noexcept { v->resize(dimensions()); }
};

template <class Var>
class NonSquareBoxDynamicDimMat : public BoxBase<NonSquareBoxDynamicDimMat<Var>, Var>,
                                  public MatBoxBase<NonSquareBoxDynamicDimMat<Var>, Var>,
                                  public NonSquareBoxCore<Var> {
 public:
  using Var_t = Var;
  using BoxBase_t = BoxBase<NonSquareBoxDynamicDimMat<Var>, Var>;
  using MatBoxBase_t = MatBoxBase<NonSquareBoxDynamicDimMat<Var>, Var>;
  static_assert(!array_traits<Var_t>::isFixedSize);
  static_assert(!array_traits<Var_t>::isVector);

  inline int dimensions() const noexcept {
    assert(this->min().size() == this->max().size());
    assert(this->min().rows() == this->max().rows());
    return this->min().size();
  }

  inline int boxRows() const noexcept {
    assert(this->min().rows() == this->max().rows());
    return this->min().rows();
  }

  inline int boxCols() const noexcept {
    assert(this->min().cols() == this->max().cols());
    return this->min().cols();
  }

  inline void setDimensions(const int _r, const int _c) noexcept {
    assert(_r > 0);
    assert(_c > 0);

    this->min().resize(_r, _c);
    this->max().resize(_r, _c);
  }

  inline void initializeSize(Var_t* v) const noexcept { v->resize(boxRows(), boxCols()); }
};

template <class Var_t>
using NonSquareBox_t = std::conditional_t<
    array_traits<Var_t>::isFixedSize,
    std::conditional_t<array_traits<Var_t>::isVector, NonSquareBoxConstDimVec<Var_t>,
                       NonSquareBoxConstDimMat<Var_t>>,

    std::conditional_t<array_traits<Var_t>::isVector, NonSquareBoxDynamicDimVec<Var_t>,
                       NonSquareBoxDynamicDimMat<Var_t>>>;

//////////////////////

}  // namespace

template <class Var, BoxShape BS = BoxShape::SQUARE_BOX>
class DiscretBox
    : public std::conditional_t<BS == BoxShape::SQUARE_BOX, SquareBox_t<Var>, NonSquareBox_t<Var>> {
 public:
  using Var_t = Var;
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;

  static_assert(!(std::is_same_v<Scalar_t, bool> && BS == BoxShape::RECTANGLE_BOX),
                "Boolean box must be a square box");
};

//////////////////////////////////

namespace {
template <class Var_t>
class SquareBoxDeltaCore {
 public:
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;
  static_assert(std::is_floating_point_v<Scalar_t>);

  inline Scalar_t& delta() noexcept { return _deltaS; }

  inline Scalar_t delta() const noexcept { return _deltaS; }

  inline Scalar_t delta(const int idx) const noexcept { return _deltaS; }

  inline void setDelta(const Scalar_t __delta) noexcept {
    assert(__delta > 0);
    _deltaS = __delta;
  }

 private:
  Scalar_t _deltaS;
};

template <class Var_t>
class NonSquareBoxDeltaCore {
 public:
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;
  static_assert(std::is_floating_point_v<Scalar_t>);

  inline Var_t& delta() noexcept { return _deltaV; }

  inline const Var_t& delta() const noexcept { return _deltaV; }

  inline Scalar_t delta(const int idx) const noexcept { return at(_deltaV, idx); }

  inline void setDelta(const Var_t& __delta) noexcept {
    for (Scalar_t deltaElement : __delta) {
      assert(deltaElement > 0);
    }
    _deltaV = __delta;
  }

 private:
  Var_t _deltaV;
};

////////////////////////////////

template <class Var_t, BoxShape BS>
class BoxContinousVec
    : public DiscretBox<Var_t, BS>,
      public std::conditional_t<BS == BoxShape::SQUARE_BOX, SquareBoxDeltaCore<Var_t>,
                                NonSquareBoxDeltaCore<Var_t>> {
 public:
  inline void setDimensions(const int dim) noexcept {
    DiscretBox<Var_t, BS>::setDimensions(dim);
    if constexpr (BS == BoxShape::RECTANGLE_BOX) {
      this->delta().resize(dim);
    }
  }
};

template <class Var_t, BoxShape BS>
class BoxContinousMat
    : public DiscretBox<Var_t, BS>,
      public std::conditional_t<BS == BoxShape::SQUARE_BOX, SquareBoxDeltaCore<Var_t>,
                                NonSquareBoxDeltaCore<Var_t>> {
 private:
  using Base = DiscretBox<Var_t, BS>;
  using conditionalBase =
      typename std::conditional_t<BS == BoxShape::SQUARE_BOX, SquareBoxDeltaCore<Var_t>,
                                  NonSquareBoxDeltaCore<Var_t>>;

 public:
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;

  inline Scalar_t delta(const int r, const int c) const noexcept {
    Base::MatBoxBase_t::assert4Size(r, c);

    return static_cast<const conditionalBase*>(this)->delta(r + c * this->boxRows());
  }

  using conditionalBase::delta;

  inline void setDimensions(const int r, const int c) noexcept {
    DiscretBox<Var_t, BS>::setDimensions(r, c);
    if constexpr (BS == BoxShape::RECTANGLE_BOX) {
      conditionalBase::delta().resize(r, c);
    }
  }
};

}  // namespace

template <class Var, BoxShape BS>
class ContinousBox : public std::conditional_t<array_traits<Var>::isVector,
                                               BoxContinousVec<Var, BS>, BoxContinousMat<Var, BS>> {
 public:
  static constexpr BoxShape Shape = BS;
  using Var_t = Var;
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;

  static_assert(std::is_floating_point_v<Scalar_t>,
                "Continous box constraint requires the type of element to be floating point types");

  inline void applyDelta(Var_t* v) const noexcept {
    assert(v->size() == this->dimensions());
    const int idx = randIdx(this->dimensions());

    at(*v, idx) += randD(-1, 1) * this->delta(idx);

    at(*v, idx) = std::min(at(*v, idx), this->max(idx));
    at(*v, idx) = std::max(at(*v, idx), this->min(idx));
  }
};

namespace {

template <class Var_t>
class GuassianBoxCore {
 public:
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;
  static_assert(std::is_floating_point_v<Scalar_t>,
                "Infinite box requries that type of element must be floating point numbers");

  inline constexpr Scalar_t min() const { return -std::numeric_limits<Scalar_t>::infinity(); }
  inline constexpr Scalar_t min(const int) const {
    return -std::numeric_limits<Scalar_t>::infinity();
  }
  inline constexpr Scalar_t max() const { return std::numeric_limits<Scalar_t>::infinity(); }
  inline constexpr Scalar_t max(const int) const {
    return std::numeric_limits<Scalar_t>::infinity();
  }

  inline Scalar_t& mean() noexcept { return _meanS; }
  inline Scalar_t mean() const noexcept { return _meanS; }

  inline Scalar_t& variance() noexcept { return _varianceS; }
  inline Scalar_t variance() const noexcept { return _varianceS; }

 private:
  Scalar_t _meanS;
  Scalar_t _varianceS;
};

}  // namespace

template <class Var>
class GaussianBox : public GuassianBoxCore<Var>, public internal::SquareBoxSizeBody<Var> {
 public:
  static constexpr BoxShape Shape = BoxShape::SQUARE_BOX;
  using Var_t = Var;
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;
  inline void initialize(Var_t* v) const noexcept {
    this->initializeSize(v);

    for (Scalar_t& val : *v) {
#warning !
    }
  }

 protected:
};

}  // namespace heu

#endif  //  HEU_BOXCONSTRAINT_HPP