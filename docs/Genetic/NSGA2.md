# NSGA2
Nondomainance sorting genetic algorithm II template class.

| Header: | `#include<HeuristicFlow/Genetic>` |
| ----: | :---- |
| Location: | [NSGA2.hpp](../../Genetic/NSGA2.hpp) |
| Inherits from: | [NSGA2Base](./NSGA2Base.md) |

<br>

## Defination
```cpp
template<typename Var_t,
         size_t ObjNum,
         DoubleVectorOption DVO,
         FitnessOption isGreaterBetter,
         RecordOption Record,
         PFOption ProtectPF,
         class ...Args> class Heu::NSGA2;
```
<br>

## Types
| Access | Name | Type | Defination |
| :----: | :----: | ----: | :---- |
| template | `DVO` | `DoubleVectorOption` |  |
| public | [`Fitness_t`](#fitness_t) | `typedef` |  |

Some other types are inherited from [NSGA2Base](./NSGA2Base.md).
<br>

## Members
All members are inherited from [NSGA2Base](./NSGA2Base.md).
<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public |  | [`NSGA2()`](#nsga2) |
| public | virtual | [`~NSGA2()`](#\~nsga2) |
<br>

## Detailed description
NSGA-II is a template for multi-objective genetic algorithm based on Pareto Optimality. It selects by nondomainance sorting, and you can order whether mutation happens on pareto front in template.

If template parameter `DVO` is `Heu::DoubleVectorOption::Eigen`, a parital specialization implementation will be activated. This speical version reloads many functions with a boosted type. **However, exposed API isn't changed.**

<br>

## Type details
### `DVO`
Defined in [Enumerations.hpp](../Global/../../Global/Enumerations.hpp):
```cpp
enum DoubleVectorOption {
    Std='S',
    Eigen='E',
    Custom='C'
};
```
This enumeration type indicates which style of double array is used as fitness value.
- `Std` refers to cpp standard arrays, including `std::valarray`(or `std::vector`) and `std::array`
- `Eigen` refers to Eigen's `Array`, including `Eigen::ArrayXd` and `Eigen::Array<double,N,1>`.
- `Custom` refers to any custom types that can behave like an array. It must support random access through `operator[]` and can get size through function `size()`.

In NSGA2, `Custom` isn't allowed. If you want to use your custom type as fitness value, use [NSGA2Base](./NSGA2Base.md).


### `Fitness_t`
Types of fitness value. It depends on the value of `ObjNum` and `DVO` as table below:

|  | `DVO==Eigen` | else |
| :----: | :----: | :----: |
| `ObjNum==Dynamic` | `Eigen::ArrayXd` | `std::valarray<double>` |
| else | `Eigen::Array<double,ObjNum,1>` | `std::array<double,ObjNum>` |

## Function details
### `NSGA2()`
Constructor.

### `~NSGA2()`
Destructor.
