
#ifndef HEU_MULTIBITSET_HPP
#define HEU_MULTIBITSET_HPP

#include <stdint.h>
#include <type_traits>
#include <memory>
#include <cmath>
#include <limits>

#include "../../Global"

//#include <iostream>

namespace heu {

template <int eleBits, typename block = size_t, class allocator_t = std::allocator<block>>
class multiBitSet {
 public:
  using block_t = block;
  static_assert(eleBits > 0);
  static_assert(eleBits <= 64);
  static_assert(std::is_integral_v<block_t>);
  static_assert(std::is_unsigned_v<block_t>);
  static_assert(sizeof(block_t) * 8 >= eleBits);

 public:
  static constexpr block_t blockMask = (block_t(1) << (eleBits)) - 1;

  static constexpr size_t blockBits = sizeof(block_t) * 8;

  using value_t = uintX_t<eleBits>;

  class reference_t;
  class iterator_t;
  class constIterator_t;

  class reference_t {
   public:
    inline operator value_t() const noexcept { return source.getValue(index); }

    inline value_t operator=(value_t newVal) noexcept { return source.setValue(index, newVal); }

    inline value_t operator=(const reference_t another) noexcept {
      return source.setValue(another.index, another);
    }

   private:
    friend class multiBitSet<eleBits>;
    friend class iterator_t;
    friend class constIterator_t;
    reference_t(multiBitSet& _src, size_t _idx) : source(_src), index(_idx) {}

    multiBitSet& source;
    size_t index;
  };

  class iterator_t {
   private:
    iterator_t(multiBitSet& src, size_t idx) : ref(src, idx) {}

    friend class multiBitSet<eleBits>;
    friend class reference_t;

    reference_t ref;

   public:
    inline reference_t& operator*() noexcept { return ref; }

    inline const reference_t& operator*() const noexcept { return ref; }

    inline iterator_t& operator++() noexcept {
      ref.index++;
      return *this;
    }

    inline iterator_t operator++(int) noexcept {
      iterator_t it;
      ref.index++;
      return it;
    }

    bool operator==(const iterator_t& another) const noexcept {
      return (&ref.source == &another.ref.source) && (ref.index == another.ref.index);
    }

    bool operator!=(const iterator_t& another) const noexcept { return !this->operator==(another); }
  };

  class constIterator_t {
   private:
    constIterator_t(multiBitSet& src, size_t idx) : ref(src, idx) {}

    friend class multiBitSet<eleBits>;
    friend class reference_t;

    reference_t ref;

   public:
    inline value_t operator*() noexcept { return ref; }

    inline value_t operator*() const noexcept { return ref; }

    inline constIterator_t& operator++() noexcept {
      ref.index++;
      return *this;
    }

    inline constIterator_t operator++(int) noexcept {
      iterator_t it;
      ref.index++;
      return it;
    }

    bool operator==(const constIterator_t& another) const noexcept {
      return (&ref.source == &another.ref.source) && (ref.index == another.ref.index);
    }
    bool operator!=(const constIterator_t& another) const noexcept {
      return !this->operator==(another);
    }
  };

 public:
  multiBitSet() {
    _data = nullptr;
    _end = nullptr;
    _size = 0;
  }

  multiBitSet(size_t __size) {
    if (__size <= 0) {
      return;
    }

    // const size_t __bytesNeed = std::ceil(float(eleBits * __size) / 8);
    const size_t __blocksNeed = std::ceil(float(__size * eleBits) / (8.0f * sizeof(block_t)));
    _data = alloc().allocate(__blocksNeed);
    _end = _data + __blocksNeed;
    _size = __size;

    memset(_data, 0, __blocksNeed);
  }

  multiBitSet(const multiBitSet& b) {
    if (b.size() > 0) {
      this->resize(b.size());
      memcpy(this->_data, b._data, b.blocks() * sizeof(block_t));
    } else {
      _data = nullptr;
      _end = nullptr;
      _size = 0;
    }
  }

  ~multiBitSet() {
    if (_data != nullptr) {
      alloc().deallocate(_data, _end - _data);
    }
  }

  inline size_t size() const noexcept { return _size; }

  inline size_t capacity() const noexcept {
    const size_t blockNum = _end - _data;
    return (blockNum * sizeof(block_t) * 8) / (eleBits);
  }

  inline block_t* data() noexcept { return _data; }

  inline size_t blocks() const noexcept { return _end - _data; }

  inline reference_t operator[](size_t idx) noexcept {
    assert(idx < size());
    return reference_t(*this, idx);
  }

  inline value_t operator[](size_t idx) const noexcept {
    assert(idx < size());
    return getValue(idx);
  }

  inline void resize(size_t newSize) noexcept {
    if (newSize <= this->capacity()) {
      _size = newSize;
      return;
    }

    const size_t newBlockNum = std::ceil(float(newSize * eleBits) / (sizeof(block_t) * 8));

    block_t* newDataPtr = alloc().allocate(newBlockNum);

    // memcpy(newDataPtr, _data, _end - _data);

    alloc().deallocate(_data, _end - _data);

    _data = newDataPtr;
    _end = _data + newBlockNum;
    _size = newSize;
  }

  inline void reserve(size_t newCap) noexcept {
    if (newCap <= this->capacity()) {
      return;
    }

    const size_t newBlockNum = std::ceil(float(newCap * eleBits) / (sizeof(block_t) * 8));

    block_t* newDataPtr = alloc().allocate(newBlockNum);

    memcpy(newDataPtr, _data, _end - _data);

    alloc().deallocate(_data, _end - _data);

    _data = newDataPtr;
    _end = _data + newBlockNum;
  }

  inline void push_back(value_t val) noexcept {
    reserve(size() + 1);
    setValue(size(), val);
    _size++;
  }

  inline void emplace_back(value_t val) noexcept { push_back(val); }

  inline void pop_back() noexcept {
    assert(_size > 0);
    _size--;
  }

  inline void clear() noexcept {
    if (_data != nullptr) {
      alloc().deallocate(_data, _end - _data);
    }

    _data = nullptr;
    _end = nullptr;
    _size = 0;
  }

  inline iterator_t begin() noexcept { return iterator_t(*this, 0); }

  inline constIterator_t begin() const noexcept { return constIterator_t(*this, 0); }

  inline iterator_t end() noexcept { return iterator_t(*this, size()); }

  inline constIterator_t end() const noexcept { return constIterator_t(*this, size()); }

 private:
  block_t* _data;
  block_t* _end;
  size_t _size;

  friend class reference_t;

  value_t setValue(const size_t index, block_t value) noexcept {
    value &= blockMask;

    const size_t firstScalarBitIdx = eleBits * index;
    const size_t lastScalarBitIdx = (index + 1) * eleBits - 1;

    block_t* blockPrevPtr = _data + firstScalarBitIdx / blockBits;
    block_t* blockNextPtr = _data + lastScalarBitIdx / blockBits;

    if (blockPrevPtr == blockNextPtr) {
      const size_t endBitIdx = ((firstScalarBitIdx / blockBits) + 1) * blockBits;
      const size_t leftMoveBits = endBitIdx - lastScalarBitIdx - 1;

      const block_t inverseMask = ~(blockMask << leftMoveBits);

      (*blockPrevPtr) &= inverseMask;

      (*blockPrevPtr) |= (value << leftMoveBits);
    } else {
      const size_t midBitIdx = (lastScalarBitIdx / blockBits) * blockBits;

      const size_t leftBitNum = midBitIdx - firstScalarBitIdx;

      const size_t rightBitNum = lastScalarBitIdx - midBitIdx + 1;

      (*blockPrevPtr) = ((*blockPrevPtr) >> leftBitNum) << leftBitNum;
      (*blockNextPtr) = ((*blockNextPtr) << rightBitNum) >> rightBitNum;

      (*blockPrevPtr) |= (value >> (rightBitNum));
      (*blockNextPtr) |= (value << (blockBits - rightBitNum));
    }

    return value;
  }

  value_t getValue(const size_t index) const noexcept {
    const size_t firstScalarBitIdx = eleBits * index;
    const size_t lastScalarBitIdx = (index + 1) * eleBits - 1;

    const block_t* blockPrevPtr = _data + firstScalarBitIdx / blockBits;
    const block_t* blockNextPtr = _data + lastScalarBitIdx / blockBits;

    if (blockPrevPtr != blockNextPtr) {
      const size_t midBitIdx = (lastScalarBitIdx / blockBits) * blockBits;

      const block_t valBeforeMask =
          ((*blockPrevPtr) << (lastScalarBitIdx - midBitIdx + 1)) |
          ((*blockNextPtr) >> (midBitIdx + blockBits - lastScalarBitIdx - 1));

      return valBeforeMask & blockMask;
    }

    else {
      const size_t endBitIdx = ((firstScalarBitIdx / blockBits) + 1) * blockBits;

      return ((*blockPrevPtr) >> (endBitIdx - lastScalarBitIdx - 1)) & blockMask;
    }
  }

  inline static allocator_t& alloc() noexcept {
    static allocator_t allo;
    return allo;
  }
};

}  // namespace heu

#endif  // HEU_MULTIBITSET_HPP