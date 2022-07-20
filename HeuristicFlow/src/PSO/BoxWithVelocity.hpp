#ifndef HEU_BOXWITHVELOCITY_HPP
#define HEU_BOXWITHVELOCITY_HPP
#include <type_traits>

#include "../../Global"

namespace heu {
namespace internal {

template <class Var_t>
class SquareBox4PSO {
 public:
  static constexpr BoxShape Shape = BoxShape::SQUARE_BOX;
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;

  inline Scalar_t& posMin() noexcept { return _posMinS; }

  inline Scalar_t posMin() const noexcept { return _posMinS; }

  inline Scalar_t& posMax() noexcept { return _posMaxS; }

  inline Scalar_t posMax() const noexcept { return _posMaxS; }

  inline void setRange(const Scalar_t min, const Scalar_t max) noexcept {
    assert(min < max);
    _posMinS = min;
    _posMaxS = max;
  }

  inline Scalar_t& maxVelocity() noexcept { return _posMaxVelocityS; }

  inline Scalar_t maxVelocity() const noexcept { return _posMaxVelocityS; }

  inline void setMaxVelocity(const Scalar_t __posMaxVelocity) noexcept {
    assert(__posMaxVelocity > 0);
    _posMaxVelocityS = __posMaxVelocity;
  }

  inline void applyConstraint4Position(Var_t* p) const noexcept {
    for (int idx = 0; idx < p->size(); idx++) {
      at(*p, idx) = std::min(_posMaxS, at(*p, idx));
      at(*p, idx) = std::max(_posMinS, at(*p, idx));
    }
  }

  inline void applyConstraint4Velocity(Var_t* v) const noexcept {
    for (int idx = 0; idx < v->size(); idx++) {
      at(*v, idx) = std::min(at(*v, idx), _posMaxVelocityS);
      at(*v, idx) = std::max(at(*v, idx), -_posMaxVelocityS);
    }
  }

 private:
  Scalar_t _posMinS;
  Scalar_t _posMaxS;
  Scalar_t _posMaxVelocityS;
};

template <class Var_t>
class NonsquareBox4PSO {
 public:
  static constexpr BoxShape Shape = BoxShape::RECTANGLE_BOX;
  using Scalar_t = typename array_traits<Var_t>::Scalar_t;

  inline Var_t& posMin() noexcept { return _posMinV; }

  inline const Var_t& posMin() const noexcept { return _posMinV; }

  inline Var_t& posMax() noexcept { return _posMaxV; }

  inline const Var_t& posMax() const noexcept { return _posMaxV; }

  inline void setRange(const Var_t& _min, const Var_t& _max) noexcept {
    if constexpr (!array_traits<Var_t>::isSizeFixed) {
      assert(_min.size() == _max.size());
    }

    if constexpr (array_traits<Var_t>::isEigenClass) {
      assert((_min.array() < _max.array()).all());
    } else {
      for (int idx = 0; idx < int(_min.size()); idx++) {
        assert(_min[idx] < _max[idx]);
      }
    }

    _posMinV = _min;
    _posMaxV = _max;
  }

  inline Var_t& maxVelocity() noexcept { return _posMaxVelocityV; }

  inline const Var_t& maxVelocity() const noexcept { return _posMaxVelocityV; }

  inline void setMaxVelocity(const Var_t& __posMaxVelocity) noexcept {
    if constexpr (array_traits<Var_t>::isEigenClass) {
      assert((__posMaxVelocity.array() > 0).all());
    } else {
      for (int idx = 0; idx < int(__posMaxVelocity.size()); idx++) {
        assert(__posMaxVelocity[idx] > 0);
      }
    }

    _posMaxVelocityV = __posMaxVelocity;
  }

  inline void applyConstraint4Position(Var_t* p) const noexcept {
    for (int idx = 0; idx < p->size(); idx++) {
      at(*p, idx) = std::min(at(_posMaxV, idx), at(*p, idx));
      at(*p, idx) = std::max(at(_posMinV, idx), at(*p, idx));
    }
  }

  inline void applyConstraint4Velocity(Var_t* v) const noexcept {
    for (int idx = 0; idx < v->size(); idx++) {
      at(*v, idx) = std::min(at(*v, idx), at(_posMaxVelocityV, idx));
      at(*v, idx) = std::max(at(*v, idx), -at(_posMaxVelocityV, idx));
    }
  }

 private:
  Var_t _posMinV;
  Var_t _posMaxV;
  Var_t _posMaxVelocityV;
};

template <class Var_t>
class SquareBoxConstDim4PSO : public SquareBox4PSO<Var_t> {
 public:
  static_assert(array_traits<Var_t>::isFixedSize);

  inline constexpr int dimensions() const noexcept { return array_traits<Var_t>::sizeCT; }
};

template <class Var_t>
class SquareBoxDynamicDim4Vec4PSO : public SquareBox4PSO<Var_t> {
 public:
  static_assert(!array_traits<Var_t>::isFixedSize);
  inline int dimensions() const noexcept { return _dimensions; }

  inline void setDimensions(const int __d) noexcept {
    assert(__d > 0);
    _dimensions = __d;
  }

 private:
  int _dimensions;
};

template <class Var_t>
class SquareBoxDynamicDim4Mat4PSO : public SquareBox4PSO<Var_t> {
 public:
  inline int boxRows() const noexcept { return _boxRows; }
  inline int boxCols() const noexcept { return _boxCols; }

  inline int dimensions() const noexcept { return _boxRows * _boxCols; }

  inline void setDimensions(const int _r, const int _c) noexcept {
    assert(_r > 0);
    assert(_c > 0);

    _boxRows = _r;
    _boxCols = _c;
  }

 private:
  int _boxRows;
  int _boxCols;
};

template <class Var_t>
using SquareBoxDynamicDim4PSO =
    typename std::conditional<array_traits<Var_t>::isMatrix, SquareBoxDynamicDim4Mat4PSO<Var_t>,
                              SquareBoxDynamicDim4Vec4PSO<Var_t>>::type;

template <class Var_t>
class NonSquareBoxConstDim4PSO : public NonsquareBox4PSO<Var_t> {
 public:
  static_assert(array_traits<Var_t>::isFixedSize);
  inline constexpr int dimensions() const noexcept { return array_traits<Var_t>::sizeCT; }
};

template <class Var_t>
class NonSquareBoxDynamicDimResizer4Vec4PSO : public NonsquareBox4PSO<Var_t> {
 public:
  inline void setDimensions(int __d) noexcept {
    assert(__d > 0);
    if constexpr (!array_traits<Var_t>::isEigenClass) {
      if (this->posMin().size() != __d) {
        this->posMin().resize(__d);
      }
      if (this->posMax().size() != __d) {
        this->posMax().resize(__d);
      }
      if (this->velocityMax().size() != __d) {
        this->velocityMax().resize(__d);
      }
    } else if constexpr (array_traits<Var_t>::isColVector) {
      if (this->posMin().size() != __d) {
        this->posMin().resize(__d, 1);
      }
      if (this->posMax().size() != __d) {
        this->posMax().resize(__d, 1);
      }
      if (this->velocityMax().size() != __d) {
        this->velocityMax().resize(__d, 1);
      }
    } else {
      if (this->posMin().size() != __d) {
        this->posMin().resize(1, __d);
      }
      if (this->posMax().size() != __d) {
        this->posMax().resize(1, __d);
      }
      if (this->velocityMax().size() != __d) {
        this->velocityMax().resize(1, __d);
      }
    }
  }
};

template <class Var_t>
class NonSquareBoxDynamicDimResizer4Mat4PSO : public NonsquareBox4PSO<Var_t> {
 private:
  inline void assert4Size() const noexcept {
    assert(this->posMin().rows() == this->posMax().rows());

    assert(this->posMax().rows() == this->velocityMax().rows());

    assert(this->posMin().cols() == this->posMax().cols());

    assert(this->posMax().cols() == this->velocityMax().cols());
  }

 public:
  static_assert(array_traits<Var_t>::isEigenClass);

  inline int boxRows() const noexcept {
    assert4Size();
    return this->posMin().rows();
  }

  inline int boxCols() const noexcept {
    assert4Size();
    return this->posMin().cols();
  }

  inline void setDimensions(const int _r, const int _c) noexcept {
    assert(_r > 0);
    assert(_c > 0);

    if (this->posMin().rows() != _r || this->posMin().cols() != _c) {
      this->posMin().resize(_r, _c);
    }

    if (this->posMax().rows() != _r || this->posMax().cols() != _c) {
      this->posMax().resize(_r, _c);
    }

    if (this->velocityMax().rows() != _r || this->velocityMax().cols() != _c) {
      this->velocityMax().resize(_r, _c);
    }
  }
};

template <class Var_t>
class NonSquareBoxDynamicDim4PSO
    : public std::conditional<array_traits<Var_t>::isVector,
                              NonSquareBoxDynamicDimResizer4Vec4PSO<Var_t>,
                              NonSquareBoxDynamicDimResizer4Mat4PSO<Var_t>>::type {
 public:
  static_assert(!array_traits<Var_t>::isFixedSize);
  inline int dimensions() const noexcept {
    assert(this->posMin().size() == this->posMax().size());
    assert(this->posMax().size() == this->velocityMax().size());

    return this->posMin().size();
  }
};

}  // namespace internal

template <class Var_t, BoxShape BS>
class Box4PSO
    : public std::conditional<
          BS == BoxShape::SQUARE_BOX,
          typename std::conditional<array_traits<Var_t>::isFixedSize,
                                    internal::SquareBoxConstDim4PSO<Var_t>,
                                    internal::SquareBoxDynamicDim4PSO<Var_t>>::type,

          typename std::conditional<array_traits<Var_t>::isFixedSize,
                                    internal::NonSquareBoxConstDim4PSO<Var_t>,
                                    internal::NonSquareBoxDynamicDim4PSO<Var_t>>::type>::type

{};

}  // namespace heu

#endif  //  HEU_BOXWITHVELOCITY_HPP