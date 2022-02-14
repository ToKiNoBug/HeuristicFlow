# OtGlobal
Empty class to put global functions and non-constant-variables.

| Header: | `#include<HeuristicFlow/Global>` |
| ----: | :---- |
| Location: | [Globals.hpp](../../Global/Globals.hpp) |

<br>

## Defination
```cpp
class OptimT::OtGlobal;
```
<br>

## Releated global variables
| Access | Type | Name | Default value |
| :----: | ----: | :---- | :----: |
| global | `std::mt19937` | [`global_mt19937`](#global_mt19937) |  |

<br>

## Related global functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| global | `inline double` | [`randD()`](#global-randd) |
| global | `inline double` | [`randD(const double min,const double max)`](#global-randdconst-double-minconst-double-max) |

<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public | `static double` | [`randD()`](#randd) |
| public | `static double` | [`randD(const double min,const double max)`](#randdconst-double-minconst-double-max) |
| public | `static uint32_t` | [`makeRandSeed()`](#makerandseed) |

<br>

## Macros
| Name | Usage |
| :----: | :---- |
| [`OT_square(x)`](#ot_squarex) | Macro function |
| [`OptimT_MAKE_GLOBAL`](#optimt_make_global) | To be expanded |

<br>

## Detailed description
This class are defined to put global functions. Besides, global variables(const and non-const) are put in the same file.

<br>

## Global Function details
### (global) `randD()`
Returns random double precision floating number in range $[0,1)$.

### (global) `randD(const double min,const double max)`
Returns random double precision floating number in range $[min,max)$.

<br>

## Member Function details
### `randD()`
**DEPRECATED**

Returns random double precision floating number in range $[0,1)$.

### `randD(const double min,const double max)`
**DEPRECATED**

Returns random double precision floating number in range $[min,max)$.

### `makeRandSeed()`
Function to provide a random seed. It's defined to provide a reliable seed. 

*previous scheme uses `std::random_device` but it's not reliable on Windows when using MingW. Thus I prefer a time-based hash value.*

If first called, it will initialize `std::srand` by a hash value of `std::time(nullptr)` and returns a seed calculated by another hash value of time.

If it's called a second time, it returns a value provided by [`global_mt19937`](#global_mt19937).

<br>

## Releated variable detailes
### `global_mt19937`
A global `std::mt19937` to be used in any condition that reqiures a `std::mt19937`, such as random shuffle or random number.

<br>

## Macro details
### `OT_square(x)`
Defination:
```cpp
#define OT_square(x) (x)*(x)
```
A macro function implementation for square.

### `OptimT_MAKE_GLOBAL`
Defination:
```cpp
#define OptimT_MAKE_GLOBAL \
std::mt19937 OptimT::global_mt19937(OptimT::OtGlobal::makeRandSeed()); \
OptimT::LogisticChaos OptimT::global_logistic(randD());
```
Packed definations of non-const global variables. **Since OptimTemplates is a header only lib, this macro must be written in a source file once. Otherwise, compile error occurs.**

Example below:

```cpp
//main.cpp
#include <OptimTemplates/Genetic>
#include <iostream>

OptimT_MAKE_GLOBAL

int main() {
    std::cout<<"random number : "<<OptimT::randD()<<std::endl;
    return 0;
}

```