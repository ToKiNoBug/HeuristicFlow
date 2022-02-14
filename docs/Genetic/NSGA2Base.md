# NSGA2Base
Direct base class of NSGA2.

| Header: | `#include<HeuristicFlow/Genetic>` |
| ----: | :---- |
| Location: | [NSGA2Base.hpp](../../Genetic/NSGA2Base.hpp) |
| Inherits from: | [NSGABase](./NSGABase.md) |
| Inherited by : | [NSGA2](./NSGA2.md) |

<br>

## Defination
```cpp
template<typename Var_t,
        size_t ObjNum,
        typename Fitness_t,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class ...Args> class Heu::NSGA2Base;
```
<br>

## Types
| Access | Name | Type | Defination |
| :----: | :----: | ----: | :---- |
| global | [`CompareOption`](#compareoption) | `enum` | `enum CompareOption : int64_t` |
| public | [`congestComposeFun`](#congestcomposefun) |`typedef` | `using congestComposeFun = double(*)(const Fitness_t *,const ArgsType*);` |
| public | [`infoUnit`](#infounit) | `struct` | `struct infoUnit` |

Some other types are inherited from [MOGABase](./MOGABase.md).
<br>

## Members
| Access | Type | Name | Default value |
| :----: | ----: | :---- | :----: |
| protected | `congestComposeFun` | [`_ccFun`](#_ccfun) | `default_ccFun_liner` |

Some other members are inherited from [MOGABase](./MOGABase.md).
<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public |  | [`NSGA2Base()`](#nsga2base) |
| public | `virtual` | [`~NSGA2Base()`](#\~nsga2base) |
| public | `void` | [`setCongestComposeFun(congestComposeFun __ccFun=default_ccFun_liner)`](#setcongestcomposefuncongestcomposefun-__ccfundefault_ccfun_liner) |
| public | `static double` | [`default_ccFun_liner(const Fitness_t * f,const ArgsType*)`](#default_ccfun_linerconst-fitness_t--fconst-argstype) |
| public | `static double` | [`default_ccFun_sphere(const Fitness_t * f,const ArgsType*)`](#default_ccfun_sphereconst-fitness_t--fconst-argstype) |
| public | `static double` | [`default_ccFun_max(const Fitness_t * f,const ArgsType*)`](#default_ccfun_maxconst-fitness_t--fconst-argstype) |
| public | `static double` | [`default_ccFun_powered<int power>(const Fitness_t * f,const ArgsType*)`](#default_ccfun_poweredint-powerconst-fitness_t--fconst-argstype) |
| public | `virtual Fitness_t` | [`bestFitness() const`](#bestfitness-const) |
| protected | `virtual void` | [`customOptWhenInitialization()`](#customoptwheninitialization) |
| protected | `static bool` | [`isStrongDomain(const Fitness_t * A,const Fitness_t * B)`](#isstrongdomainconst-fitness_t--aconst-fitness_t--b) |
| protected | `static bool` | [`universialCompareFun<int64_t>(const infoUnit *,const infoUnit *)`](#universialcomparefunint64_tconst-infounit-const-infounit-) |
| protected | `virtual void` | [`select()`](#select) |

Some other functions are inherited from [MOGABase](./MOGABase.md).
<br>

## Macros
| Name | Usage |
| :----: | :----: |
| [`Heu_NSGA2_DO_PARALLELIZE`](#Heu_nsga2_do_parallelize) | conditional compiling |

<br>

## Detailed description
This class prodives a fundamental implementation of NSGA2 but it's not so easy-usable, and it doesn't supports vectorized boosting.

<br>

## Type details
### `CompareOption`
Defination:
```cpp
enum CompareOption : int64_t {
    CompareByCongestion=-1,
    CompareByDominantedBy=-2
};
```

### `congestComposeFun`
Function pointer to calculate congestion.

### `infoUnit`
```cpp
struct Heu::NSGA2Base::infoUnit
    : public Heu::NSGABase::infoUnitBase
{
public:
    Fitness_t congestion;
};
```
This struct is used to store nondomainance sorting-related informations in select operation. A pointer to this struct, `infoUnit*`, has access to every information to a gene. 

There's several sorting by according to different attribute, thus a vector of `infoUnit*` is widely used when sorting.

`congestion` stores the partial congestion in each objective. That's why it has the same type with fitness but it's not a fitness value. 
   
   **Mention that when composing congestion, the first element of `congestion` is occupied to store the total congestion and other elements aren't used.**

<br>

## Member details
### `_ccFun`
Function pointer to compose congestion value. 

After partial congestion on each objective is calculated, it's required to compose then into a congestion value(`double`).

<br>

## Function details
### `NSGA2Base()`
Default constructor.

### `~NSGA2Base()`
Default destructor.

### `setCongestComposeFun(congestComposeFun __ccFun=default_ccFun_liner)`
Set the value of [`_ccFun`](#_ccfun).

It's not a must to call this function when using NSGA2, instead, it's an optional one.

### `default_ccFun_liner(const Fitness_t * f,const ArgsType*)`
A candidate function for [`_ccFun`](#_ccfun). It composes congestion by adding difference on each objective linearly. It's the default function.
$$
{Congestion} = \sum _{i=1}^{ObjNum}{\Delta f_i}
$$
To some extent, it's like **Manhattan distance**.

### `default_ccFun_sphere(const Fitness_t * f,const ArgsType*)`
A candidate function for [`_ccFun`](#_ccfun). It compoeses congestion like a sphere using following formula:
$$
{Congestion} = \sqrt{\sum _{i=1}^{ObjNum}{\Delta f_i}^2}
$$
To some extent, it's like **Eculidean distance**.


### `default_ccFun_max(const Fitness_t * f,const ArgsType*)`
A candidate function for [`_ccFun`](#_ccfun). It select the maximun parital congestion among all objectives as the total congestion.

Formula:
$$
{Congestion} = \max _{1\leq i\leq ObjNum} {\Delta f_i}
$$
To some extent, it's like **Chebyshev distance**.

### `default_ccFun_powered<int power>(const Fitness_t * f,const ArgsType*)`
A candidate function for [`_ccFun`](#_ccfun).

Formula:
$$
{Congestion} = {(\sum _{i=1}^{ObjNum}{\Delta f_i}^p)}^{\frac{1}{p}} \\
\text{Where }p=\text{power in the template.}
$$
To some extent, it's like **Minkowski distance**.


### `bestFitness() const`
Implementation for pure virtual function defined in [GABase](./GABase.md). It's defined to provide a `Fitness_t` to be recorded. This function will compose the best value of each objective. 

**Usually such a fitness value can domainate the whole population, however corresponding solution rarely exists.**

If you aren't satisfied to such recording method, inherit `NSGA2` and do it yourself.

### `customOptWhenInitialization()`
Reimplementation to this function defined in [GABase](./GABase.md). 

It clears `_pfGenes`, tells it to reserve space for twice the size of population size. Besides, it assigns [`prevFrontSize`](#prevfrontsize) to -1 and congestion composing function to [`default_ccFun_liner(const Fitness_t * f,const ArgsType*)`](#default_ccfun_linerconst-fitness_t--fconst-argstype).

### `isStrongDomain(const Fitness_t * A,const Fitness_t * B)`
Returns whether A strong domainates B.

### `universialCompareFun<int64_t>(const infoUnit *,const infoUnit *)`
Defination:
```cpp
template<int64_t objIdx>
    static bool universialCompareFun(const infoUnit * A,const infoUnit * B);
```
Template parameter `objIdx` is used to indicate how to compare 2 `const infoUnit *`.

If `objIdx` equals to `CompareByCongestion`(-1), it's a comparison function to sort according to total congestion value in descending order. If some genes in a pareto layer will be selected and others in this layer will be eliminated, genes with greater congestion value has greater chance to be selected.

If `objIdx` equals to `CompareByDominatedBy`(-2), it's a comparison function to sort according to `infoUnit::domainedByNum` in an acsending order. This sorting step will seperate the whole population into several pareto layers.

If `objIdx` is an non-negative integer, it's a comparison function to sort according to the `infoUnit::sortIdx`-th objective value. It doesn't matter whether such sorting is in an acsending order or a decsending order.

Some template metaprogramming tricks will be used to generate an static constant array of function pointers consists of different instances of this function template with parameter ranging from 0 to `ObjNum-1`.

If objective number is dynamic, `ObjNum-1` will be replaced by 255(see [`Heu_MOGA_RTObjNum_MaxObjNum`](MOGABase.md)). 

I believe in most cases no body will tries to solve a multi-objective problem with more than 255 objectives using NSGA-II.

### `select()`
This is the core of NSGA2. It applies non-dominated sorting in following steps:
1. Make a `std::vector` of `infoUnit` named `pop` that each element corresponds to an element in population(`std::list<Gene>`).
2. Make a `std::vector` of `infoUnit*` named `sortSpace` that each pointer element points to an element in `pop`.
3. Calculate `infoUnit::domainedByNum` for `pop`.
4. Sort elements in `sortSpace` according to `infoUnit::domainedByNum`.
5. Make a `std::list<std::vector<infoUnit *>>` named `paretoLayers` and insert each element in `sortSpace`. The whole population is devided into several layers. For better performance, each `std::vector<infoUnit *>` will reserve space for the whole population. Here we use vector since sorting might happen in a single layer.
6. Make a `std::queue` of `infoUnit *` named `selected` to store all selected genes. Insert pareto front and each layers into `selected` until a layer can't be completedly inserted(I name this layer *K*). Congestion are used to select several genes in *K*.
7. Sort the whole population according to every objective values to calculate partial congestion on each objective. Each objective will be sorted once.
8. Compose partial congestion into total congestion for layer *K*.
9. Sort layer *K* according to parital congestion.
10. Inserts each element in *K* until size of `selected` equals to user assigned population size(`GABase::_option.populationSize`).
11. Mark every selected genes and erase unmarked genes from population.

Every sorting steps use `std::sort` and different comparision function are used to sort according to different attributes of a gene.

Really complicated, isn't it?


<br>

## Macro details
### `Heu_NSGA2_DO_PARALLELIZE`
This macro should only be defined if you defined `Heu_DO_PARALLELIZE` before including HeuristicFlow and you want some multi-threading boostings for `NSGA2`'s selection operator.