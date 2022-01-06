# GABase
Abstrcat base class of genetic algorithm solvers.

| Header: | `#include<Genetic>` |
| ----: | :---- |
| Location: | [GABase.h](./../../GA/GABase.h) |
| Inherited by : | [SOGA](./SOGA.md) |

<br>

## Defination
```cpp
template<typename Var_t,typename Fitness_t,bool Record,class ...Args> class OptimT::GABase;
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
| protected | `std::vector<Fitness_t>` | [`_record`](#_record) |  |
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
| public | `const std::vector<Fitness_t>` | [`record() const`](#record-const) |
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
GABase implemented basic initialization, fitness calculation, crossover, mutation and run functions. However, since all operators in GA has lots relationship with encoding method, these operators must be implemented by users, and user should pass these implementation via function pointers.




<br>

## Type details
### `Var_t`
Encoded solution type. For instance, when solving TSP problems, what we really want to solve is a permulation to arrange order between nodes, but the permulation is encoded into a double array and decoded via sorting. So in this instance, `Var_t` is assigned to `std::vector<double>` instead of an permulation type. You should decode when calculating fitness value.


### `Fitness_t`
Besides, since type of fitness value can be any type, which made **comparsion between unknown fitness type is impossible** in template, **`select()` function is pure virtual and GABase is abstract**. A derived class SOGA inherits GABase and implemented `select()` function with `double` as fitness value type. So `SOGA` is a single-object genetic algorithm solver, that's why I name it SOGA.

By the way, this also leaves room for multi-object optimizators.

### `Record`
This boolean template value determines whether instance record the changes of fitness value or not. If true, you call an extra partial special version of GABase with an extra member [`_record`](#_record) and function [`record() const`](#record-const).

### `Args...`
Extra custom parameters. They can be some datas to calculate the fitness of a Var_tiable like trainning data set, or some other parameters like maximum and minimum value of `Var_t`, or even trainning rate. Also you can leave it blank. Several examples:
```cpp
OptimT::GABase<std::array<double,16>,pathL_t,true,DotsData_t> * tspSolver;
OptimT::GABase<BPNetwork_t,loss_t,false,DataSet_t> * networkTrainner;
```

### `ArgsType`
`Args...` will be packed in a tuple (`std::tuple`) and such tuple are called **args**, which means parameters of this optimization task.

### `initializeFun`
Function pointer type to initialize an individual in population. This function will be called only when initializing. For best flexibility, a constant pointer to parameters is provided.

### `fitnessFun`
Fucntion pointer type to calculate fitness for any individual.

This function will be called to each individual in each generation, so it often spend most time when optimizing. When you are solving a complex problem, make this function as fast as possible.

### `crossoverFun`
Function pointer type to do crossover. The first 2 constant pointers can be called *parent1* and *parent2* while call rest non-constant pointers *child1* and *child2*. When doing crossover, parents are chosen from population while children are newly inserted to population. This function should assign children's value by parents' value.

Note that children are constructed with **default constructor**, **make sure that your `Var_t` type does have a default constructor**.

### `mutateFun`
Function pointer type to do mutation. It has same type with `initializeFun` but they have different meanings. This function should slightly modify given `Var_t`. Same as above, these functions have constant access to parameters.

### `Gene`
Gene type encapsulating a `Var_t` object. It has following members (and functions):
| Type/Return Type | Name | Type |
| ----: | :---- | :----: |
| `Var_t` | `self` | public member |
| `bool` | `_isCalculated` | public member |
| `double` | `_Fitness` | public member |
| `void` | `setUncalculated()` | public function |
| `bool`  | `isCalculated() const` | public function |
| `double` | `fitness() const` | public function |

To avoid repeated fitness caculation, an extra boolean is added to note whether a `Gene` has already get its fitness value. You can call `void Gene::setUncalculated()` to set a Gene to uncalculated.

### `otherOptFun`
It means "other operate function". This function is called in every generation, after selection but before crossover. In this function, you can modify not only args but also the whole population. In most condition this step is not needed, but I love flexibility.
   
**All these function pointers can be lambda. When constructed, all these function pointers has an empty lambda as default value.**

### `GeneIt`
It's same as iterator type of stl list containing `Gene` types.

<br>

## Member details
### `_eliteIt`
### `_population`
### `_option`
### `_generation`
### `_failTimes`
### `_record`
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
### `record() const`
### `option() const`
### `generation() const`
### `failTimes() const`
### `args() const`
### `calculateAll()`
### `select()=0`
### `crossover()`
### `mutate()`

