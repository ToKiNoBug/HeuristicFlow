#include <algorithm>
#include "../../Global"

namespace heu {

namespace {
template <class Derived, class Var_t>
class BoxBase {
  // using Var_t = typename Derived::Var_t;
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;

 public:
  inline void initialize(Var_t* v) const noexcept {
    static_cast<Derived*>(this)->initializeSize(v);
    for (int idx = 0; idx < v->size(); idx++) {
      if constexpr (std::is_same_v<Scalar_t, bool>) {
        at(*v, idx) = randIdx(2) % 2 == 0;
      } else if constexpr (std::is_integral_v<Scalar_t>) {
        at(*v, idx) = randIdx<Scalar_t>(static_cast<Derived*>(this)->min(idx),
                                        1 + static_cast<Derived*>(this)->max(idx));
      } else {
        at(*v, idx) = randD(static_cast<Derived*>(this)->min(idx),
                            static_cast<const Derived*>(this)->max(idx));
      }
    }
  }

  inline void applyConstraint(Var_t* v) const noexcept {
    for (int idx = 0; idx < v->size(); idx++) {
      at(*v, idx) = std::min(at(*v, idx), static_cast<const Derived*>(this)->max(idx));
      at(*v, idx) = std::max(at(*v, idx), static_cast<const Derived*>(this)->min(idx));
    }
  }

 protected:
  inline void assert4Size(const int idx) noexcept {
    assert(idx >= 0 && idx < static_cast<Derived*>(this)->dimensions());
  }
};

////////////////

template <class Derived, class Var_t>
class MatBoxBase {
 public:
 protected:
  inline void assert4Size(const int r, const int c) noexcept {
    assert(r >= 0 && r < static_cast<Derived*>(this)->boxRows());
    assert(c >= 0 && c < static_cast<Derived*>(this)->boxCols());
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

////////////////////////////

template <class Var>
class SquareBoxConstDimVec : public BoxBase<SquareBoxConstDimVec<Var>, Var>,
                             public SquareBoxCore<Var> {
 public:
  using Var_t = Var;
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;
  static_assert(array_traits<Var_t>::isFixedSize);
  static_assert(array_traits<Var_t>::isVector);

  inline constexpr int dimensions() const noexcept { return array_traits<Var_t>::sizeCT; }

  inline void initializeSize(Var_t*) const noexcept {}

 private:
};

template <class Var>
class SquareBoxConstDimMat : public BoxBase<SquareBoxConstDimMat<Var>, Var>,
                             public MatBoxBase<SquareBoxConstDimMat<Var>, Var>,
                             public SquareBoxCore<Var> {
 public:
  using Var_t = Var;
  static_assert(array_traits<Var_t>::isFixedSize);
  static_assert(!array_traits<Var_t>::isVector);

 public:
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;

  inline constexpr int dimensions() const noexcept { return array_traits<Var_t>::sizeCT; }

  inline constexpr int boxRows() const noexcept { return array_traits<Var_t>::rowsCT; }

  inline constexpr int boxCols() const noexcept { return array_traits<Var_t>::colsCT; }

  inline void initializeSize(Var_t*) const noexcept {}

 private:
};

template <class Var>
class SquareBoxDynamicDimVec : public BoxBase<SquareBoxDynamicDimVec<Var>, Var>,
                               public SquareBoxCore<Var> {
 public:
  using Var_t = Var;
  static_assert(!array_traits<Var_t>::isFixedSize);
  static_assert(array_traits<Var_t>::isVector);

  inline int dimensions() const noexcept { return _dimensions; }

  inline void setDimensions(const int dim) noexcept {
    assert(dim > 0);
    _dimensions = dim;
  }

 private:
  int _dimensions;
};

template <class Var>
class SquareBoxDynamicDimMat : public BoxBase<SquareBoxDynamicDimMat<Var>, Var>,
                               public MatBoxBase<SquareBoxDynamicDimMat<Var>, Var>,
                               public SquareBoxCore<Var> {
 public:
  using Var_t = Var;
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
  static_assert(array_traits<Var_t>::isFixedSize);
  static_assert(array_traits<Var_t>::isVector);

  inline constexpr int dimensions() const noexcept { return array_traits<Var_t>::sizeCT; }
};

template <class Var>
class NonSquareBoxConstDimMat : public BoxBase<NonSquareBoxConstDimMat<Var>, Var>,
                                public MatBoxBase<NonSquareBoxConstDimMat<Var>, Var>,
                                public NonSquareBoxCore<Var> {
 public:
  using Var_t = Var;

  static_assert(array_traits<Var_t>::isFixedSize);
  static_assert(!array_traits<Var_t>::isVector);

  inline constexpr int dimensions() const noexcept { return array_traits<Var_t>::sizeCT; }

  inline constexpr int boxRows() const noexcept { return array_traits<Var_t>::rowsCT; }

  inline constexpr int boxCols() const noexcept { return array_traits<Var_t>::colsCT; }
};

template <class Var>
class NonSquareBoxDynamicDimVec : public BoxBase<NonSquareBoxDynamicDimVec<Var>, Var>,
                                  public NonSquareBoxCore<Var> {
 public:
  using Var_t = Var;

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
};

template <class Var>
class NonSquareBoxDynamicDimMat : public BoxBase<NonSquareBoxDynamicDimMat<Var>, Var>,
                                  public MatBoxBase<NonSquareBoxDynamicDimMat<Var>, Var>,
                                  public NonSquareBoxCore<Var> {
 public:
  using Var_t = Var;

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

template <class Var_t, BoxShape BS = BoxShape::SQUARE_BOX>
using DiscretBox =
    std::conditional_t<BS == BoxShape::SQUARE_BOX, SquareBox_t<Var_t>, NonSquareBox_t<Var_t>>;

}  // namespace heu