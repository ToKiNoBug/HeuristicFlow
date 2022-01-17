# NSGA2
Nondomainance sorting genetic algorithm II template class.

| Header: | `#include<Genetic>` |
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
| global | [`PFOption`](#pfoption) | `enum` | `PFOption : unsigned char{...}` |
| template | [`ObjNum`](#objnum) | `size_t` |  |
| public | [`This_t`](#this_t) | `typedef` | `using This_t = NSGA2;` |
| public | [`Base_t`](#base_t) |`typedef` | `using Base_t = GABase<Var_t,std::array<double,ObjNum>,Record,Args...>;` |
| public | [`Fitness_t`](#fitness_t) | `typedef` | `using Fitness_t = std::array<double,ObjNum>;` |
| public | [`congestComposeFun`](#congestcomposefun) |`typedef` | `using congestComposeFun = double(*)(const Fitness_t *,const ArgsType*);` |
| public | [`infoUnit`](#infounit) | `struct` | `struct infoUnit` |

Some other types are inherited from [GABase](./GABase.md).
<br>

## Members
| Access | Type | Name | Default value |
| :----: | ----: | :---- | :----: |
| protecte | `size_t` | [`prevFrontSize`](#prevfrontsize) |  |
| protecte | `std::unordered_set<const Gene*>` | [`_pfGenes`](#_pfgenes) |  |
| protecte | `congestComposeFun` | [`_ccFun`](#_ccfun) | `default_ccFun_liner` |

Some other members are inherited from [GABase](./GABase.md).
<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public |  | [`NSGA2()`](#nsga2) |
| public |  | [`~NSGA2()`](#\~nsga2) |
| public | `void` | [`setCongestComposeFun(congestComposeFun __ccFun=default_ccFun_liner)`](#setcongestcomposefuncongestcomposefun-__ccfundefault_ccfun_liner) |
| public | `static double` | [`default_ccFun_liner(const Fitness_t * f,const ArgsType*)`](#default_ccfun_linerconst-fitness_t--fconst-argstype) |
| public | `static double` | [`default_ccFun_sphere(const Fitness_t * f,const ArgsType*)`](#default_ccfun_sphereconst-fitness_t--fconst-argstype) |
| public | `static double` | [`default_ccFun_max(const Fitness_t * f,const ArgsType*)`](#default_ccfun_maxconst-fitness_t--fconst-argstype) |
| public | `static double` | [`default_ccFun_powered<int power>(const Fitness_t * f,const ArgsType*)`](#default_ccfun_poweredint-powerconst-fitness_t--fconst-argstype) |
| public | `virtual Fitness_t` | [`bestFitness() const`](#bestfitness-const) |
| public | `void` | [`paretoFront(std::vector<Fitness_t> & front)`](#paretofrontstdvectorfitness_t--front) |
| public | `void` | [`paretoFront(std::vector<std::pair<const Var_t*,const Fitness_t*>> & front)`](#paretofrontstdvectorstdpairconst-var_tconst-fitness_t--front) |
| public | `const std::unordered_set<const Gene*> &` | [`pfGenes() const`](#pfgenes-const) |
| protected | `virtual void` | [`customOptWhenInitialization()`](#customoptwheninitialization) |
| protected | `static bool` | [`isStrongDomain(const Fitness_t * A,const Fitness_t * B)`](#isstrongdomainconst-fitness_t--aconst-fitness_t--b) |
| protected | `static bool` | [`compareFun_DomainedBy(const infoUnit * A,const infoUnit * B)`](#comparefun_domainedbyconst-infounit--aconst-infounit--b) |
| protected | `static bool` | [`compareFun_Fitness(const infoUnit * A,const infoUnit * B)`](#comparefun_fitnessconst-infounit--aconst-infounit--b) |
| protected | `static bool` | [`compareFun_Congestion(const infoUnit * A,const infoUnit * B)`](#comparefun_congestionconst-infounit--aconst-infounit--b) |
| protected | `virtual void` | [`select()`](#select) |
| protected | `virtual void` | [`mutate()`](#mutate) |

Some other functions are inherited from [GABase](./GABase.md).
<br>

## Detailed description
NSGA-II is a template for multi-objective genetic algorithm based on Pareto Optimality. It selects by nondomainance sorting, and you can order whether mutation happens on pareto front in template.

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
Number of objectives. If a integer less than 2 provided, static assertion emerged.

### `This_t`
An shortcut to type of NSGA2 class. Equals to `typeof(*this)`.

### `Base_t`
An shortcut to base class.

### `Fitness_t`
An shortcut to `std::array<double,ObjNum>`.

### `congestComposeFun`
Function pointer to calculate congestion.

### `infoUnit`
```cpp
struct OptimT::NSGA2::infoUnit
{
    public:
        bool isSelected;
        size_t sortIdx;
        size_t domainedByNum;
        GeneIt_t iterator;
        Fitness_t congestion;
};
```
This struct is used to store nondomainance sorting-related informations in select operation. A pointer to this struct, `infoUnit*`, has access to every information to a gene. 

<br>

## Member details
### `prevFrontSize`
Stores the size of pareto front in previous generation.

### `_pfGenes`
A unordered set to mark the pareto front.

### `_ccFun`
Function pointer to compose congestion value. 

After partial congestion on each objective is calculated, it's required to compose then into a congestion value(`double`).

<br>

## Function details
### `NSGA2()`
Constructor.

### `~NSGA2()`
Destructor.

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

### `paretoFront(std::vector<Fitness_t> & front)`
Get pareto front consisted of a vector of `Fitness_t`.

### `paretoFront(std::vector<std::pair<const Var_t*,const Fitness_t*>> & front)`
Get pareto front consisted of a vector of `<std::pair<const Var_t*,const Fitness_t*>`. It provided both solution(`Var_t`) and objective value(`Fitness_t`).

### `pfGenes() const`
Returns a const reference to member [`_pfGenes`](#_pfgenes).

### `customOptWhenInitialization()`
### `isStrongDomain(const Fitness_t * A,const Fitness_t * B)`
### `compareFun_DomainedBy(const infoUnit * A,const infoUnit * B)`
### `compareFun_Fitness(const infoUnit * A,const infoUnit * B)`
### `compareFun_Congestion(const infoUnit * A,const infoUnit * B)`
### `select()`
### `mutate()`


