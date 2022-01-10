# MatrixFixedSize
Col-major fixed-size 2d array template implementation.

| Header: | `#include<SimpleMatrix>` |
| ----: | :---- |
| Location: | [MatrixFixedSize.hpp](../../SimpleMatrix/MatrixFixedSize.hpp) |

<br>

## Defination
```cpp
template<class Scalar_t,size_t Rows,size_t Cols> class MatrixFixedSize;
```
<br>

## Types
| Access | Name | Type | Defination |
| :----: | :----: | ----: | :---- |
|  | `Scalar_t` | template `typename` |  |
|  | `Rows` | template `size_t` |  |
|  | `Cols` | template `size_t` |  |
| public | `iterator` | `typedef` | `using iterator = Scalar_t*;` |
| public | `citerator` |`typedef` | `using citerator = const Scalar_t*;` |

<br>

## Members
| Access | Type | Name | Default value |
| :----: | ----: | :---- | :----: |
| protected | `Scalar_t` | [`array[Rows*Cols];`](#array) |  |

<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public |  | `MatrixFixedSize()` |
| public |  | `~MatrixFixedSize()` |
| public | `explicit` | `MatrixFixedSize(const MatrixFixedSize & src)` |
| public | `iterator` | `begin()` |
| public | `iterator` | `end()` |
| public | `static size_t` | `size()` |
| public | `static size_t` | `rows()` |
| public | `static size_t` | `cols()` |
| public | `Scalar_t &` | `operator()(size_t n)` |
| public | `Scalar_t &` | `operator()(size_t r,size_t c)` |
| public | `Scalar_t *` | `data()` |
| public | `const Scalar_t *` | `cdata()` |
| public | `const MatrixFixedSize &` | `operator=(const MatrixFixedSize & src)` |

<br>

## Detailed description
A simple fixed-size 2d array of any type, without math computation. It is a col-major one.

<br>

## Type details
### `Scalar_t`
Type of element.

### `Rows`
Number of row;

### `Cols`
Number of coloumn.

### `iterator`
STL-like Iterator.

### `citerator`
STL-like constant iterator.

<br>

## Member details
### `array`
A 2d array of `Scalar_t` to store elements.

<br>

## Function details
### `MatrixFixedSize()`
Default constructor.

### `~MatrixFixedSize()`
Default destructor.

### `MatrixFixedSize(const MatrixFixedSize & src)`
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

