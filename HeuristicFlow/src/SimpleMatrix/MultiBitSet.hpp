
#ifndef HEU_MULTIBITSET_HPP
#define HEU_MULTIBITSET_HPP

#include <stdint.h>
#include <type_traits>
#include <memory>
#include <cmath>
#include <limits>

#include "../../Global"

#include <iostream>

namespace heu {

/**
 * \ingroup HEU_SIMPLEMATRIX
 * \brief A shrinked array of N-bit unsigned integers.
 *
 * \tparam eleBits Bits of each element
 * \tparam block A unsinged integer type to store the data, default value is `size_t`
 * \tparam allocator_t Type of memory allocator. Default value is `std::allocator`
 */
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
  /// Type of element
  using value_t = uintX_t<eleBits>;
  /// The maximum value of each element
  static constexpr block_t blockMask = (block_t(1) << (eleBits)) - 1;
  /// Bits of block
  static constexpr size_t blockBits = sizeof(block_t) * 8;

  class reference_t;
  class iterator_t;
  class constIterator_t;

  /**
   * \brief Type of reference. Since N-bit unsigned integer is not a fundamental type, pointer and
   * native reference to such type is invalid. `operator value()` and `operator=` are reloaded to
   * enable reading and writing.
   *
   */
  class reference_t {
   public:
    /**
     * \brief Convert to a value.
     *
     * \return value_t The value of that element.
     */
    inline operator value_t() const noexcept { return source->getValue(index); }

    /**
     * \brief Assign new value to the referenced element. If new value is greater than the maximum,
     * extra bits will be ignored.
     *
     * \param newVal The new value
     * \return value_t Type of value
     */
    inline value_t operator=(value_t newVal) noexcept { return source->setValue(index, newVal); }

    /**
     * \brief Assign a new value with another reference.
     *
     * \param another Another reference
     * \return value_t Type of value.
     */
    inline value_t operator=(const reference_t another) noexcept {
      return source->setValue(index, value_t(another));
    }

    /**
     * \brief If the reference is invalid
     *
     * \return true If the reference is valid
     * \return false If the reference is invalid
     */
    inline bool isNull() const noexcept {
      return (source == nullptr) || (index < 0) || (index >= source->size());
    }

    /*
        template <class T>
        inline value_t operator+=(T newVal) noexcept {
          return operator=(operator value_t() + newVal);
        }

        template <class T>
        inline value_t operator-=(T newVal) noexcept {
          return operator=(operator value_t() - newVal);
        }

        template <class T>
        inline value_t operator*=(T newVal) noexcept {
          return operator=(operator value_t() * newVal);
        }

        template <class T>
        inline value_t operator/=(T newVal) noexcept {
          return operator=(operator value_t() / newVal);
        }
        */

   private:
    friend class multiBitSet<eleBits>;
    friend class iterator_t;
    friend class constIterator_t;
    reference_t() : source(nullptr), index(0){};
    reference_t(multiBitSet* _src, size_t _idx) : source(_src), index(_idx) {}

    inline void swap(reference_t& another) noexcept {
      std::swap(this->source, another.source);
      std::swap(this->index, another.index);
    }

    inline void evaluateFrom(const reference_t& another) noexcept {
      this->source = another.source;
      this->index = another.index;
    }

    multiBitSet* source;
    size_t index;
  };

  template <class it_t>
  class iteratorBase {
   public:
    /**
     * \brief Iterate to the next element and return the iterator AFTER moving.
     * Called like ++it
     *
     * \return const it_t& The iterator after moving.
     */
    inline const it_t& operator++() noexcept {
      static_cast<it_t*>(this)->ref.index++;
      return static_cast<const it_t&>(*this);
    }

    /**
     * \brief Iterate to the next element and return the element BEFORE moving.
     * Called like it++
     *
     * \return it_t The iterator before moving
     */
    inline it_t operator++(int) noexcept {
      it_t it = static_cast<it_t&>(*this);
      static_cast<it_t*>(this)->ref.index++;
      return it;
    }

    /**
     * \brief Iterate to the previous element and return the element AFTER moving.
     * Called like --it
     *
     * \return const it_t& The iterator after moving
     */
    inline const it_t& operator--() noexcept {
      static_cast<it_t*>(this)->ref.index--;
      return static_cast<const it_t&>(*this);
    }

    /**
     * \brief Iterate to the previous element and return the element BEFORE moving.
     * Called like it--
     *
     * \return it_t The iterator before moving
     */
    inline it_t operator--(int) noexcept {
      it_t it = static_cast<it_t&>(*this);
      static_cast<it_t*>(this)->ref.index--;
      return it;
    }

    inline it_t operator+(int64_t idx) const noexcept {
      return it_t(static_cast<const it_t*>(this)->ref.source,
                  static_cast<const it_t*>(this)->ref.index + idx);
    }

    inline it_t operator-(int64_t idx) const noexcept {
      return it_t(static_cast<const it_t*>(this)->ref.source,
                  static_cast<const it_t*>(this)->ref.index - idx);
    }

    inline it_t operator[](int64_t idx) const noexcept { return *this + idx; }

    inline size_t index() const noexcept { return static_cast<const it_t*>(this)->ref.index; }

    inline bool isNull() const noexcept { return static_cast<const it_t*>(this)->ref.isNull(); }

   protected:
    inline bool isSame(const reference_t& anotherRef) const noexcept {
      return (static_cast<const it_t*>(this)->ref.source == anotherRef.source) &&
             (static_cast<const it_t*>(this)->ref.index == anotherRef.index);
    }
  };

  /**
   * \brief Type of iterator. `operator*, `operator++(int)`, `operator--`, `operator--(int)`
   * and `operator++`, have been reloaded.
   *
   */
  class iterator_t : public iteratorBase<iterator_t> {
   private:
    iterator_t(multiBitSet* src, size_t idx) : ref(src, idx) {}

    friend class multiBitSet<eleBits>;
    friend class reference_t;
    friend class iteratorBase<iterator_t>;

    reference_t ref;

   public:
    /**
     * \brief Convert an iterator to a reference. (acutally a native reference of reference_t)
     *
     * \return reference_t& A reference to the pointed element.
     */
    inline reference_t& operator*() noexcept { return ref; }

    /**
     * \brief Convert an iterator to a read only reference.
     *
     * \return const reference_t& A read only reference to the pointed element.
     */
    inline const reference_t& operator*() const noexcept { return ref; }

    /**
     * \brief Compare between different iterators
     *
     * \param another Another iterator
     * \return true When iterators are same
     * \return false When iterators are different
     */
    template <class anotherIt_t>
    inline bool operator==(const anotherIt_t& another) const noexcept {
      return this->isSame(another.ref);
    }

    /**
     * \brief Compare between different iterators
     *
     * \param another Another iterator
     * \return true When iterators are different
     * \return false When iterators are same.
     */
    template <class anotherIt_t>
    inline bool operator!=(const anotherIt_t& another) const noexcept {
      return !(this->isSame(another.ref));
    }

    /**
     * \brief Copy the value of iterator
     *
     */
    const iterator_t& operator=(const iterator_t& another) noexcept {
      this->ref.index = another.ref.index;
      this->ref.source = another.ref.source;
      return *this;
    }
  };

  /**
   * \brief The read only iterator. `operator*, `operator++(int)`, `operator--`, `operator--(int)`
   * and `operator++`, have been reloaded.
   *
   */
  class constIterator_t : public iteratorBase<constIterator_t> {
   private:
    constIterator_t(multiBitSet* src, size_t idx) : ref(src, idx) {}

    friend class multiBitSet<eleBits>;
    friend class reference_t;
    friend class iteratorBase<constIterator_t>;

    reference_t ref;

   public:
    /**
     * \brief Construct a new constIterator_t object from non-const iterator
     *
     */
    constIterator_t(const iterator_t& parent) {
      this->ref.index = parent.ref.index;
      this->ref.source = parent.ref.source;
    }

    /**
     * \brief Get the value of pointed element.
     *
     * \return value_t The value of pointed element.
     */
    inline value_t operator*() const noexcept { return ref; }

    /**
     * \brief Compare between different iterators
     *
     * \param another Another iterator
     * \return true When iterators are same
     * \return false When iterators are different
     */
    template <class anotherIt_t>
    inline bool operator==(const anotherIt_t& another) const noexcept {
      return this->isSame(another.ref);
    }

    /**
     * \brief Compare between different iterators
     *
     * \param another Another iterator
     * \return true When iterators are different
     * \return false When iterators are same.
     */
    template <class anotherIt_t>
    inline bool operator!=(const anotherIt_t& another) const noexcept {
      return !(this->isSame(another.ref));
    }

    /**
     * \brief Copy the const iterator
     *
     */
    const constIterator_t& operator=(const constIterator_t& another) noexcept {
      this->ref.index = another.ref.index;
      this->ref.source = another.ref.source;
      return *this;
    }
  };

 public:
  /**
   * \brief Default initializer
   *
   */
  multiBitSet() {
    _data = nullptr;
    _end = nullptr;
    _size = 0;
  }

  /**
   * \brief Initialize with given size
   *
   * \param __size The initial size.
   */
  multiBitSet(size_t __size) {
    if (__size <= 0) {
      return;
    }

    // const size_t __bytesNeed = std::ceil(float(eleBits * __size) / 8);
    const size_t __blocksNeed = std::ceil(float(__size * eleBits) / (blockBits));
    _data = alloc().allocate(__blocksNeed);
    _end = _data + __blocksNeed;
    _size = __size;

    memset(_data, 0, __blocksNeed * sizeof(block_t));
  }

  /**
   * \brief Copy constructor (deep copy)
   *
   * \param b
   */
  multiBitSet(const multiBitSet& b) {
    _data = nullptr;
    _end = nullptr;
    _size = 0;
    if (b.size() > 0) {
      this->resize(b.size());

      memcpy(this->_data, b._data, b.blocks() * sizeof(block_t));
    }
  }

  /**
   * \brief Move constructor (shallow copy)
   *
   */
  multiBitSet(multiBitSet&& b) {
    _data = b._data;

    _end = b._end;
    _size = b._size;
  }

  /**
   * \brief Destroy the multi Bit Set object
   *
   */
  ~multiBitSet() {
    if (_data != nullptr) {
      alloc().deallocate(_data, _end - _data);
    }
  }

  /**
   * \brief The element amount
   *
   * \return size_t The amount of elements.
   */
  inline size_t size() const noexcept { return _size; }

  /**
   * \brief Get the capacity
   *
   * \return size_t The amount of elements that multiBitSet can contain without allocating memory
   */
  inline size_t capacity() const noexcept {
    const size_t blockNum = _end - _data;
    return (blockNum * sizeof(block_t) * 8) / (eleBits);
  }

  /**
   * \brief The pointer of data
   *
   * \return block_t* A pointer to data
   */
  inline block_t* data() noexcept { return _data; }

  /**
   * \brief Blocks uesd to store all elements.
   *
   * \return size_t Number of blocks
   */
  inline size_t blocks() const noexcept { return _end - _data; }

  inline reference_t at(size_t idx) noexcept { return this->operator[](idx); }

  inline value_t at(size_t idx) const noexcept { return this->operator[](idx); }

  inline reference_t operator[](size_t idx) noexcept {
    assert(idx < size());
    return reference_t(this, idx);
  }

  inline value_t operator[](size_t idx) const noexcept {
    assert(idx < size());
    return getValue(idx);
  }

  /**
   * \brief Change the size of multiBitSet
   *
   * \param newSize
   */
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

  /**
   * \brief Reserve memory for elements
   *
   * \param newCap
   */
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

  inline void shrink_to_fit() noexcept {
    const size_t blocksNeed = std::ceil(float(size() * eleBits) / blockBits);

    if (blocksNeed == blocks()) {
      return;
    }

    block_t* newData = alloc().allocate(blocksNeed);

    memcpy(newData, _data, sizeof(block_t) * blocksNeed);

    alloc().deallocate(_data, blocks());

    _data = newData;
    _end = _data + blocksNeed;
  }

  /**
   * \brief Add a element to the back of array.
   */
  inline void push_back(value_t val) noexcept {
    reserve(size() + 1);
    setValue(size(), val);
    _size++;
  }

  /**
   * \brief Add a element to the back of array.
   */
  inline void emplace_back(value_t val) noexcept { push_back(val); }

  /**
   * \brief Popout he last element.
   *
   */
  inline void pop_back() noexcept {
    assert(_size > 0);
    _size--;
  }

  /**
   * \brief Clear all elements.
   *
   */
  inline void clear() noexcept {
    if (_data != nullptr) {
      alloc().deallocate(_data, _end - _data);
    }

    _data = nullptr;
    _end = nullptr;
    _size = 0;
  }

  /**
   * \brief Deep copy
   *
   * \param another Copy source
   * \return const multiBitSet& A reference to self
   */
  inline const multiBitSet& operator=(const multiBitSet& another) noexcept {
    if (another.size() <= 0) {
      resize(0);
      return *this;
    }

    resize(another.size());
    const size_t blocksNeed = std::ceil(float(size() * eleBits) / blockBits);

    memcpy(_data, another._data, blocksNeed * sizeof(block_t));
    return *this;
  }

  inline reference_t front() noexcept { return reference_t(this, 0); }

  inline value_t front() const noexcept { return getValue(0); }

  inline reference_t back() noexcept { return reference_t(this, size() - 1); }

  inline value_t back() const noexcept { return getValue(size() - 1); }

  inline iterator_t begin() noexcept { return iterator_t(this, 0); }

  inline constIterator_t begin() const noexcept { return constIterator_t(this, 0); }

  inline iterator_t end() noexcept { return iterator_t(this, size()); }

  inline constIterator_t end() const noexcept { return constIterator_t(this, size()); }

 private:
  block_t* _data;
  block_t* _end;
  size_t _size;

  friend class reference_t;

  value_t setValue(const size_t index, block_t value) noexcept {
    value &= blockMask;  // remove extra bits

    const size_t firstScalarBitIdx = eleBits * index;
    const size_t lastScalarBitIdx = (index + 1) * eleBits - 1;
    // a element may exists in 2 blocks. Compute addresses of both blocks
    block_t* blockPrevPtr = _data + firstScalarBitIdx / blockBits;
    block_t* blockNextPtr = _data + lastScalarBitIdx / blockBits;

    // if pointers are the same, the element is in a single block
    if (blockPrevPtr == blockNextPtr) {
      const size_t endBitIdx = ((firstScalarBitIdx / blockBits) + 1) * blockBits;
      const size_t leftMoveBits = endBitIdx - lastScalarBitIdx - 1;

      const block_t inverseMask = ~(blockMask << leftMoveBits);

      (*blockPrevPtr) &= inverseMask;

      (*blockPrevPtr) |= (value << leftMoveBits);
    } else {
      // the element is stored in two blocks

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