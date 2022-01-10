# OtGlobal
Empty class to put global functions and non-constant-variables.

| Header: | `#include<OtGlobal>` |
| ----: | :---- |
| Location: | [Globals.hpp](../../Global/Globals.hpp) |

<br>

## Defination
```cpp
class OptimT::OtGlobal;
```
<br>

## Releated variables
| Access | Type | Name | Default value |
| :----: | ----: | :---- | :----: |
| global | `std::random_device` | [`global_random_device`](#global_random_device) |  |
| global | `std::mt19937` | [`global_mt19937`](#global_mt19937) |  |

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public | `static double` | [`randD()`](#randd) |
| public | `static double` | [`randD(const double min,const double max)`](#randdconst-double-minconst-double-max) |

<br>

## Macros
| Name | Usage |
| :----: | :---- |
| [`OT_square(x)`](#ot_squarex) | Macro function |
| [`OPTIMT_MAKE_GLOBAL`](#optimt_make_global) | To be expanded |

<br>

## Detailed description
This class are defined to put global functions. Besides, global variables(const and non-const) are put in the same file.

<br>

## Function details
### `randD()`
Returns random double precision floating number in range $[0,1)$.

### `randD(const double min,const double max)`
Returns random double precision floating number in range $[min,max)$.

<br>

## Releated variable detailes
### `global_random_device`
A global `std::random_device` to be used in any condition that reqiures a `std::random_device`, such as random shuffle or random number.

### `global_mt19937`
A global `std::mt19937` to be used in any condition that reqiures a `std::mt19937`, such as random shuffle or random number.

## Macro details
### `OT_square(x)`
Defination:
```cpp
#define OT_square(x) (x)*(x)
```
A macro function implementation for square.

### `OPTIMT_MAKE_GLOBAL`
Defination:
```cpp
#define OPTIMT_MAKE_GLOBAL \
std::random_device OptimT::global_random_device; \
std::mt19937 OptimT::global_mt19937(OptimT::global_random_device());
```
Packed definations of non-const global variables. **Since OptimTemplates is a header only lib, this macro must be written in a source file once. Otherwise, compile error occurs.**

Example below:

```cpp
//main.cpp
#include <OptimTemplates/Genetic>
#include <iostream>

OPTIMT_MAKE_GLOBAL

int main() {
    std::cout<<"random number : "<<OptimT::OtGlobal::randD()<<std::endl;
    return 0;
}

```