#ifndef HEU_POOLVECTOR_HPP
#define HEU_POOLVECTOR_HPP

#include <type_traits>
#include <stdint.h>
#include <memory>

#include "InternalHeaderCheck.h"

namespace heu {

template <class Scalar_t, class allocator_t = std::allocator<Scalar_t>>
class poolVector {
 private:
  Scalar_t* _data;
  size_t _size;
  size_t _capacity;

  static_assert(!std::is_same<Scalar_t, void>::value, "Scalar_t mustn't be void!");
  static const bool isClass = std::is_class<Scalar_t>::value;

  inline static allocator_t& alloc() {
    static allocator_t alloctor;
    return alloctor;
  }

 public:
  poolVector() {
    _data = nullptr;
    _size = 0;
    _capacity = 0;
  }

  poolVector(const poolVector&) = delete;
  poolVector(poolVector&&) = delete;

  poolVector(size_t sz) {
    _data = nullptr;
    _size = 0;
    _capacity = 0;
    resize(sz);
  }

  ~poolVector() {
    if (_data != nullptr) {
      if (isClass)
        for (size_t i = 0; i < _capacity; i++) {
          alloc().destroy(_data + i);
        }
      alloc().deallocate(_data, _capacity);
    }
  }

  inline Scalar_t* data() noexcept { return _data; }
  inline const Scalar_t* data() const noexcept { return _data; }
  inline size_t size() const noexcept { return _size; }
  inline size_t capacity() const noexcept { return _capacity; }

  inline Scalar_t& operator[](size_t idx) {
    assert(idx < _size);
    return _data[idx];
  }
  inline const Scalar_t& operator[](size_t idx) const {
    assert(idx < _size);
    return _data[idx];
  }

  inline Scalar_t& operator()(size_t idx) {
    assert(idx < _size);
    return _data[idx];
  }
  inline const Scalar_t& operator()(size_t idx) const {
    assert(idx < _size);
    return _data[idx];
  }

  inline Scalar_t& at(size_t idx) {
    assert(idx < _size);
    return _data[idx];
  }
  inline const Scalar_t& at(size_t idx) const {
    assert(idx < _size);
    return _data[idx];
  }

  Scalar_t* begin() noexcept { return _data; }
  Scalar_t* end() noexcept { return _data + _size; };
  const Scalar_t* begin() const noexcept { return _data; }
  const Scalar_t* end() const noexcept { return _data + _size; };

  Scalar_t& front() {
    assert(_size > 0);
    return _data[0];
  }
  const Scalar_t& front() const {
    assert(_size > 0);
    return _data[0];
  }
  Scalar_t& back() {
    assert(_size > 0);
    return _data[_size - 1];
  }
  const Scalar_t& back() const {
    assert(_size > 0);
    return _data[_size - 1];
  }

  inline void reserve(size_t newCapacity) {
    if (newCapacity <= capacity()) {
      return;
    }

    Scalar_t* newPtr = alloc().allocate(newCapacity);
    if (isClass) {
      for (size_t idx = 0; idx < newCapacity; idx++) {
        alloc().construct(newPtr + idx);
      }
    }

    for (size_t idx = 0; idx < _size; idx++) {
      newPtr[idx] = _data[idx];
    }

    if (isClass) {
      for (size_t idx = 0; idx < _capacity; idx++) {
        alloc().destroy(_data + idx);
      }
    }

    alloc().deallocate(_data, _capacity);

    _data = newPtr;
    _capacity = newCapacity;
  }

  inline void resize(size_t newSize) {
    if (newSize != size()) {
      if (newSize > _capacity) {
        if (_data != nullptr) {
          if (isClass)
            for (size_t i = 0; i < _capacity; i++) {
              alloc().destroy(_data + i);
            }
          alloc().deallocate(_data, _capacity);
        }

        _data = alloc().allocate(newSize);
        _capacity = newSize;

        if (isClass)
          for (size_t i = 0; i < _capacity; i++) {
            alloc().construct(_data + i);
          }
      }
    }
    _size = newSize;
  }

  inline void clear() { _size = 0; }

  inline void pop_back() {
    assert(_size > 0);
    _size--;
  }

  inline void push_back(const Scalar_t& back) {
    if (_size < _capacity) {
      _data[_size] = back;
      _size++;
    } else {
      const size_t newSize = _size + 1;
      const size_t newCapacity = newSize;
      Scalar_t* newPtr = alloc().allocate(newCapacity);
      // construct new objects
      if (isClass) {
        for (size_t idx = 0; idx < newCapacity; idx++) {
          alloc().construct(newPtr + idx);
        }
      }
      // copy old objects to new place
      for (size_t idx = 0; idx < _size; idx++) {
        newPtr[idx] = _data[idx];
      }
      // copy the last object
      newPtr[newSize - 1] = back;
      // destroy the old objects
      if (isClass) {
        for (size_t idx = 0; idx < _capacity; idx++) {
          alloc().destroy(_data + idx);
        }
      }

      alloc().deallocate(_data, _capacity);
      _data = newPtr;
      _size = newSize;
      _capacity = newCapacity;
    }
  }

  // inline void emplace_back(const Scalar_t& back) { emplace_back(std::move(back)); }
};

}  //  namespace heu

#endif  //  HEU_POOLVECTOR_HPP