# GABase
Abstrcat base class of genetic algorithm solvers.

| Header: | `#include<Genetic>` |
| ----: | :---- |
| Location: | [GABase.h](./../../GA/GABase.h) |
| Inherited by : | [SOGA](./SOGA.md) |

<br>

## Defination
```cpp
template<typename Var_t,typename Fitness_t,class ...Args> class OptimT::GABase;
```
<br>

## Types
| Access | Name | Type | Defination |
| :----: | :----: | ----: | :---- |
|  | [`Var_t`](#var_t) | `template` |  |
|  | [`Fitness_t`](#fitness_t) | `template` |  |
|  | [`Args...`](#args...) | `template` |  |
| public | [`ArgsType`](#argstype) | `typedef` | `typedef std::tuple<Args...> ArgsType;` |
| public | [`initializeFun`](#initializefun) | `typedef` |  `typedef void(* initializeFun)(Var_t*,const ArgsType*);` |
| public | [`fitnessFun`](#fitnessfun) | `typedef` | `typedef Fitness_t(* fitnessFun)(const Var_t*,const ArgsType*);` |
| public | [`crossoverFun`](#crossoverfun) | `typedef` | `typedef void(* crossoverFun)(const Var_t*,const Var_t*,Var_t*,Var_t*,const ArgsType*);` |
| public | [`mutateFun`](#mutatefun) | `typedef` | `typedef initializeFun mutateFun;` |
| public | [`Gene`](#gene) | `class` | `class Gene` |
| public | [`otherOptFun`](#otheroptfun) | `typedef` | `typedef void(* otherOptFun)(ArgsType*,std::list<Gene>*,size_t generation,size_t failTimes,const GAOption*);` |
| protected | [`GeneIt`](#geneit) | `typedef` | `typedef typename std::list<Gene>::iterator GeneIt;` |

<br>

## Members
| Access | Type | Name | Default value |
| :----: | ----: | :---- | :----: |
| protected | `GeneIt` | [`_eliteIt`](#_eliteit) |  |
| protected | `std::list<Gene>` | [`_population`](#_population) |  |
| protected | [`GAOption`](./GAOption.md) | [`_option`](#_option) |  |
| protected | `size_t` | [`_generation`](#_generation) |  |
| protected | `size_t` | [`_failTimes`](#_failtimes) |  |
| protected | `std::vector<Fitness_t>` | [`_recording`](#_recording) |  |
| protected | `ArgsType` | [`_args`](#_args) |  |
| protected | `initializeFun` | [`_initializeFun`](#_initializefun) | empty lambda |
| protected | `fitnessFun` | [`_fitnessFun`](#_fitnessfun) | empty lambda |
| protected | `crossoverFun` | [`_crossoverFun`](#_crossoverfun) | empty lambda |
| protected | `mutateFun` | [`_mutateFun`](#_mutatefun) | empty lambda |
| protected | `otherOptFun` | [`_otherOptFun`](#_otheroptfun) | empty lambda |

<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public |  | [`GABase()`](#gabase) |
| public | `virtual` | [`~GABase()`](#\~gabase) |
| public | `virtual void` | [`run()`](#run) |
| public | `std::list<Gene>::iterator` | [`eliteIt() const`](#eliteit-const) |
| public | `const Var_t &` | [`result() const`](#result-const) |
| public | `const std::list<Gene> &` | [`population() const`](#population-const) |
| public | `const std::vector<Fitness_t>` | [`recording() const`](#recording-const) |
| public | `const GAOption &` | [`option() const`](#option-const) |
| public | `size_t generation()` | [`generation() const`](#generation-const) |
| public | `size_t failTimes()` | [`failTimes() const`](#failtimes-const) |
| public | `const ArgsType &` | [`args() const`](#args-const) |
| protected | `virtual void` | [`calculateAll()`](#calculateall) |
| protected | `virtual void` | [`select()=0`](#select0) |
| protected | `virtual void` | [`crossover()`](#crossover) |
| protected | `virtual void` | [`mutate()`](#mutate) |

<br>

## Detailed description
Class/Strcut description here.

<br>

## Type details
### `Var_t`
### `Fitness_t`
### `Args...`
### `ArgsType`
### `initializeFun`
### `fitnessFun`
### `crossoverFun`
### `mutateFun`
### `Gene`
### `otherOptFun`
### `GeneIt`

<br>

## Member details
### `_eliteIt`
### `_population`
### `_option`
### `_generation`
### `_failTimes`
### `_recording`
### `_args`
### `_initializeFun`
### `_fitnessFun`
### `_crossoverFun`
### `_mutateFun`
### `_otherOptFun`

<br>

## Function details
### `GABase()`
### `~GABase()`
### `run()`
### `eliteIt() const`
### `result() const`
### `population() const`
### `recording() const`
### `option() const`
### `generation() const`
### `failTimes() const`
### `args() const`
### `calculateAll()`
### `select()=0`
### `crossover()`
### `mutate()`

