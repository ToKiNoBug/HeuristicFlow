# MOGAAbstract
An abstract base class for multiple objective genetic algorithm solvers.

| Header: | `#include<OptimTemplates/Genetic>` |
| ----: | :---- |
| Location: | [MOGAAbstract.hpp](../../Genetic/MOGAAbstract.hpp) |
| Inherits from: | [GABase](GABase.md) |
| Inherited by : | MOGABase,NSGA2Base,[NSGA2](NSGA2.md) |

<br>

## Defination
```cpp
template<typename Var_t,
        size_t ObjNum,
        typename Fitness_t,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class ...Args> class MOGAAbstract;
```
<br>

## Types
| Access | Name | Type | Defination |
| :----: | :----: | ----: | :---- |
| global | [`PFOption`](#pfoption) | `enum` | `PFOption : unsigned char{...}` |
| template | [`ObjNum`](#objnum) | `size_t` |  |
| template | [`pfOpt`](#pfopt) | `PFOption` |  |
| public | [`Base_t`](#base_t) |`typedef` | `using Base_t = GABase<Var_t,Fitness_t,Record,Args...>;` |

Some other types are inherited from [GABase](./GABase.md).
<br>

## Members
| Access | Type | Name | Default value |
| :----: | ----: | :---- | :----: |
| protected | `size_t` | [`prevFrontSize`](#prevfrontsize) |  |
| protected | `size_t` | [`prevPFCheckSum`](#prevPFCheckSum) |  |
| protected | `std::unordered_set<const Gene*>` | [`_pfGenes`](#_pfgenes) |  |

Some other members are inherited from [GABase](./GABase.md).
<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public |  | [`MOGAAbstract()`](#mogaabstract) |
| public | `virtual` | [`~MOGAAbstract()`](#\~mogaabstract) |
| public | `void` | [`paretoFront(std::vector<Fitness_t> & front)`](#paretofrontstdvectorfitness_t--front) |
| public | `void` | [`paretoFront(std::vector<std::pair<const Var_t*,const Fitness_t*>> & front)`](#paretofrontstdvectorstdpairconst-var_tconst-fitness_t--front) |
| public | `const std::unordered_set<const Gene*> &` | [`pfGenes() const`](#pfgenes-const) |
| protected | `virtual void` | [`mutate()`](#mutate) |

Some other functions are inherited from [GABase](./GABase.md).
<br>


## Detailed description
It's a base class that contained functions to implement a Pareto front. A hash calculating virtual function is implemented which enables derived classess to compare whether the Pareto front is changed or not in a generation.

<br>

## Type details
### `PFOption`
```cpp
enum PFOption : unsigned char {
    PARETO_FRONT_DONT_MUTATE=true,
    PARETO_FRONT_CAN_MUTATE=false
};
```
This enumeration type indicates that whether the pareto front will be protected from mutation when algorithm is running.

### `ObjNum`
Number of objectives. If a integer equals to 1 provided, static assertion emerged. Use `OptimT::Dynamic`(`0`) to mark that number of objectives are determined at runtime.

### `pfOpt`
Whether to protect the Pareto front from mutation or not.

### `Base_t`
An shortcut to direct base class.


<br>

## Member details
### `prevFrontSize`
Stores the size of pareto front in previous generation.

### `prevPFCheckSum`
Hash checksum of Pareto front in previous generation.

### `_pfGenes`
A unordered set to mark the pareto front.


<br>

## Function details
### `MOGAbstract()`
Default constructor.

### `~MOGAbstract()`
Default destrcutor.

### `paretoFront(std::vector<Fitness_t> & front)`
Get Pareto front consisted of a vector of `Fitness_t`.

### `paretoFront(std::vector<std::pair<const Var_t*,const Fitness_t*>> & front)`
Get Pareto front consists of a vector of `std::pair<const Var_t*,const Fitness_t*>`. It provides both solution(`Var_t`) and objective value(`Fitness_t`).

### `pfGenes() const`
Returns a const reference to member [`_pfGenes`](#_pfgenes).

### `mutate()`
Implementation of pure virtual function defined in [GABase](./GABase.md). It's slightly different from mutation operation in [SOGA](./SOGA.md). 

If template parameter [`PFOption`](#pfoption) is assigned to `PARETO_FRONT_DONT_MUTATE`(`true`), the pareto front (elements in [`_pfGenes`](#_pfgenes)) won't mutate. Otherwise every genes in population have probability to mutate.
