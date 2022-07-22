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

#ifndef HEU_SIZEBODY4BOXCONSTRAINT_HPP
#define HEU_SIZEBODY4BOXCONSTRAINT_HPP

#include "../../Global"

namespace heu {

namespace {
template <class Var_t>
class SquareBoxSizeBody4ConstDimVec {
 public:
  inline void initializeSize(Var_t*) const noexcept {}

  inline constexpr int dimensions() const noexcept { return array_traits<Var_t>::sizeCT; }
};

template <class Var_t>
class SquareBoxSizeBody4ConstDimMat : public SquareBoxSizeBody4ConstDimVec<Var_t> {
 public:
  inline constexpr int boxRows() const noexcept { return array_traits<Var_t>::rowsCT; }

  inline constexpr int boxCols() const noexcept { return array_traits<Var_t>::colsCT; }
};

template <class Var_t>
class SquareBoxSizeBody4DynamicDimVec {
 public:
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
  int _dimensions;
};

template <class Var_t>
class SquareBoxSizeBody4DynamicDimMat {
 public:
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

}  // namespace

namespace internal {

template <class Var>
using SquareBoxSizeBody = std::conditional_t<
    array_traits<Var>::isVector,
    std::conditional_t<array_traits<Var>::isFixedSize, SquareBoxSizeBody4ConstDimVec<Var>,
                       SquareBoxSizeBody4DynamicDimVec<Var>>,
    std::conditional_t<array_traits<Var>::isFixedSize, SquareBoxSizeBody4ConstDimMat<Var>,
                       SquareBoxSizeBody4DynamicDimMat<Var>>>;

}  // namespace internal

}  // namespace heu

#endif  //  HEU_SIZEBODY4BOXCONSTRAINT_HPP