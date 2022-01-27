# MOGABase
Base class for multiple objective genetic solvers.

| Header: | `#include<OptimTemplates/Genetic>` |
| ----: | :---- |
| Location: | [MOGABase.hpp](../../Genetic/MOGABase.hpp) |
| Inherits from: | [MOGAAbstract](./MOGAAbstract.md) |
| Inherited by : | NSGA2Base, [NSGA2](./NSGA2.md) |

<br>

## Defination
```cpp
template<typename Var_t,
        size_t ObjNum,
        typename Fitness_t,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class ...Args> class MOGABase;
```
<br>

## Types
All types are inherited from [MOGAAbstract](./MOGAAbstract.md).
<br>

## Members
| Access | Type | Name | Default value |
| :----: | ----: | :---- | :----: |
| protected | `size_t` | [`_objectiveNum`](#_objectivenum) |  |

Some other members are inherited from [MOGAAbstract](./MOGAAbstract.md).
<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public |  | [`MOGABase()`](#mogabase) |
| public | `virtual` | [`~MOGABase()`](#\~mogabase) |
| public | `inline size_t` | [`objectiveNum() const`](#objectivenum-const) |
| public | `void` | [`setObjectiveNum(size_t)`](#setobjectivenumsize_t) |

<br>

## Detailed description
This layer of inheritance aims to support solvers with fixed and dynamic objective numbers. This template class has a default implementation for fixed objective numbers and a partial specialization for dynamic.

Lots of code-copying can be avoided by such inheritance and parital specialization.

<br>

## Member details
### `_objectiveNum`
Numbers of objectives. 

**This member only exists when template parameter `ObjNum` is `OptimT::Dynamic`(`0`).**

<br>

## Function details
### `MOGABase()`
Default constructor.

### `~MOGABase()`
Default destrcutor.

### `objectiveNum() const`
Get number of objectives. If template parameter `ObjNum` is `OptimT::Dynamic`(`0`), this function returns the value of member [`_objectiveNum`](#_objectivenum), otherwise it returns the value of `ObjNum`.

### `setObjectiveNum(size_t)`
Set the number of objectives.

**This function only exists when template parameter `ObjNum` is `OptimT::Dynamic`(`0`).**