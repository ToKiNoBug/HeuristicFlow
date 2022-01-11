# MatrixDynamicSize
Col-major dynamic-size 2d array template implementation.

| Header: | `#include<SimpleMatrix>` |
| ----: | :---- |
| Location: | [MatrixDynamicSize.hpp](../../SimpleMatrix/MatrixDynamicSize.hpp) |

<br>

## Defination
```cpp
template<class Scalar_t> class OptimT::MatrixDynamicSize;
```
<br>

## Types
| Access | Name | Type | Defination |
| :----: | :----: | ----: | :---- |
|  | [`Scalar_t`](#scalar_t) | template `typename` |  |
| public | [`iterator`](#iterator) | `typedef` | `using iterator = Scalar_t*;` |
| public | [`citerator`](#citerator) |`typedef` | `using citerator = const Scalar_t*;` |

<br>

## Members
| Access | Type | Name | Default value |
| :----: | ----: | :---- | :----: |
| protected | `Scalar_t *` | [`dataPtr`](#dataptr) | `nullptr` |
| protected | `size_t` | [`rowNum`](#rownum) | 0 |
| protected | `size_t *` | [`colNum`](#colnum) | 0 |

<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public |  | [`MatrixDynamicSize()`](#matrixfixedsize) |
| public |  | [`~MatrixDynamicSize()`](#\~matrixfixedsize) |
| public | `explicit` | [`MatrixDynamicSize(const MatrixDynamicSize &)`](#matrixdynamicsizeconst-matrixdynamicsize--src) |
| public | `iterator` | [`begin()`](#begin) |
| public | `iterator` | [`end()`](#end) |
| public | `size_t` | [`size()`](#size) |
| public | `size_t` | [`rows()`](#rows) |
| public | `size_t` | [`cols()`](#cols) |
| public | `void` | [`resize(size_t r,size_t c)`](#resizesize_t-rsize_t-c) |
| public | `Scalar_t &` | [`operator()(size_t n)`](#operatorsize_t-n) |
| public | `Scalar_t &` | [`operator()(size_t r,size_t c)`](#operatorsize_t-rsize_t-c) |
| public | `Scalar_t *` | [`data()`](#data) |
| public | `const Scalar_t *` | [`cdata()`](#cdata) |
| public | `const MatrixDynamicSize &` | [`operator=(const MatrixDynamicSize &)`](#operatorconst-matrixfixedsize--src) |

<br>

## Detailed description
A simple dynamic-size 2d array of any type, without math computation. It is a col-major one.

<br>

## Type details
### `Scalar_t`
Type of element.

### `iterator`
STL-like Iterator.

### `citerator`
STL-like constant iterator.

<br>

## Member details
### `dataPtr`
Data pointer to the first element.

### `rowNum`
Number of rows.

### `colNum`
Number of coloums.

<br>

## Function details
### `MatrixFixedSize()`
Default constructor. Matrix size will be set to zero so no memory allocation in this function.

### `~MatrixFixedSize()`
Default destructor. Memory will be freed.

### `MatrixDynamicSize(const MatrixDynamicSize & src)`
Deep-copy constructor.

### `begin()`
STL-like begin function that return an iterator to the first element.

### `end()`
STL-like end function that return an iterator after the last element.

### `size()`
Count of elements.

### `rows()`
Count of rows.

### `cols()`
Count of coloumns.

### `resize(size_t r,size_t c)`
Resize matrix. If element count is not changed, it beheaves like conservative resize; otherwise every element will be set to undefined value.

### `operator()(size_t n)`
Random access to the n-th element.

### `operator()(size_t r,size_t c)`
Random access.

### `data()`
Writeable data pointer.

### `cdata()`
Readonly data pointer.

### `operator=(const MatrixFixedSize & src)`
Deep copy function.

