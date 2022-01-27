# NSGA2
Nondomainance sorting genetic algorithm II template class.

| Header: | `#include<OptimTemplates/Genetic>` |
| ----: | :---- |
| Location: | [NSGA2.hpp](../../Genetic/NSGA2.hpp) |
| Inherits from: | [GABase](./GABase.md) |

<br>

## Defination
```cpp
template<typename Var_t,size_t ObjNum,
         FitnessOption isGreaterBetter,
         RecordOption Record,
         PFOption ProtectPF,
         class ...Args> class OptimT::NSGA2;
```
<br>

## Types
| Access | Name | Type | Defination |
| :----: | :----: | ----: | :---- |
| public | [`Fitness_t`](#fitness_t) | `typedef` | `using Fitness_t = std::array<double,ObjNum>;` |

Some other types are inherited from [GABase](./GABase.md).
<br>

## Members
All members are inherited from [GABase](./GABase.md).
<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public |  | [`NSGA2()`](#nsga2) |
| public | virtual | [`~NSGA2()`](#\~nsga2) |
<br>

## Detailed description
NSGA-II is a template for multi-objective genetic algorithm based on Pareto Optimality. It selects by nondomainance sorting, and you can order whether mutation happens on pareto front in template.


<br>

## Type details
### `Base_t`
An shortcut to base class.

### `Fitness_t`
An shortcut to `std::array<double,ObjNum>`.

## Member details

<br>

## Function details
### `NSGA2()`
Constructor.

### `~NSGA2()`
Destructor.
