# NSGABase
Base class for nondominated sorting series of genetic algorithm solvers, like NSGA2 and NSGA3.

| Header: | `#include<OptimTemplates/Genetic>` |
| ----: | :---- |
| Location: | [NSGABase.hpp](../../Genetic/NSGABase.hpp) |
| Inherits from: | [MOGABase](./MOGABase.md) |
| Inherited by : | [NSGA2Base](./NSGA2Base.md), [NSGA2](./NSGA2.md) |

<br>

## Defination
```cpp
template<typename Var_t,
        size_t ObjNum,
        typename Fitness_t,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class ...Args> class OptimT::NSGABase;
```
<br>

## Types
| Access | Name | Type | Defination |
| :----: | :----: | ----: | :---- |
| public | [`infoUnitBase`](#infounitbase) | `class` | `infoUnitBase` |

Some other types are inherited from [MOGABase](./MOGABase.md).
<br>

## Members
All members are inherited from [MOGABase](./MOGABase.md).

<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public |  | [`NSGABase()`](#nsgabase) |
| public | `virtual` | [`~NSGABase()`](#\~nsgabase) |
| protected | `virtual void` | [`calculateDominatedNum(infoUnitBase ** pop,const size_t popSizeBefore)`](#calculatedominatednuminfounitbase--popconst-size_t-popsizebefore) |
| protected | `void` | [`updatePF(const infoUnitBase ** pfs,const size_t curFrontSize)`](#updatepfconst-infounitbase--pfsconst-size_t-curfrontsize) |

Some other functions are inherited from [MOGABase](./MOGABase.md).
<br>

## Macros
| Name | Usage |
| :----: | :----: |
| [`OptimT_MAKE_NSGABASE_TYPES`](#optimt_make_nsgabase_types) | to be expanded |

<br>

## Detailed description
Base class for nondominated sorting series of genetic algorithm solvers, like NSGA2 and NSGA3. This class implemented some basical function used in NS.

<br>

## Type details
### `infoUnitBase`
Defination:
```cpp
struct OptimT::NSGABase::infoUnitBase
{
public:
    bool isSelected;
    size_t domainedByNum;
    GeneIt_t iterator;
};
```
This struct stores basic informations of a gene used in nondominated sorting.

1. `isSelected` marks whether a gene is selected to form the next generation.
2. `domainedByNum` means the number of genes in previous that strong domainates this gene. **Value 0 means it's a member of the pareto front.**
3. `iterator` is a `std::list<Gene>::iterator` to a gene.

<br>

## Function details
### `NSGABase()`
Default constructor.

### `~NSGABase()`
Default destructor.

### `calculateDominatedNum(infoUnitBase ** pop,const size_t popSizeBefore)`
Function to calculate `infoUnitBase::domainedByNum` of a population. `pop` is an array of `infoUnit*` and `popSizeBefore` is the population size before selection.

### `updatePF(const infoUnitBase ** pfs,const size_t curFrontSize)`
Function to update Pareto front. If element count in PF doesn't change, checksum of PF will be computed to compare if PF has changed.

<br>

## Macro details
### `OptimT_MAKE_NSGABASE_TYPES`
Defination:
```cpp
#define OptimT_MAKE_NSGABASE_TYPES \
OptimT_MAKE_GABASE_TYPES \
using infoUnitBase_t = typename Base_t::infoUnitBase;
```
Replacement for `OptimT_MAKE_GABASE_TYPES` used in derived class of this.