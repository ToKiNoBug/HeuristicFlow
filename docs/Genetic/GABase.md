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
| protected | [`GAOption`](./GAOption.md) | [`_option`] |  |
| protected | `size_t` | `_generation` |  |
| protected | `size_t` | `_failTimes` |  |
| protected | `std::vector<Fitness_t>` | `_recording` |  |
| protected | `ArgsType` | `_args` |  |
| protected | `fitnessFun` | `_fitnessFun` | empty lambda |
| protected | `initializeFun` | `_initializeFun` | empty lambda |
| protected | `crossoverFun` | `_crossoverFun` | empty lambda |
| protected | `mutateFun` | `_mutateFun` | empty lambda |
| protected | `otherOptFun` | `_otherOptFun` | empty lambda |

<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public |  | `GAOption()` |

<br>

## Detailed description
Class/Strcut description here.

<br>

## Type details
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
### 
### 
### 
### 
### 
### 
### 
### 
### 
### 
### 
### 

<br>

## Function details
### function name
### function name

