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
| protected | [`GeneIt_t`](#GeneIt_t) | `typedef` | `typedef typename std::list<Gene>::iterator GeneIt_t;` |

<br>

## Members
| Access | Type | Name | Default value |
| :----: | ----: | :---- | :----: |
| protected | `GeneIt_t` | [`_eliteIt`](#_eliteit) |  |
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
| public | `virtual void` | [`initialize(initializeFun,fitnessFun,crossoverFun,mutateFun,otherOptFun=nullptr,const GAOption&=GAOption(),const ArgsType&=ArgsType())`](#initialize) |
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

Use `RECORD_FITNESS`(`true`) here to enable recording function, or use `DONT_RECORD_FITNESS`(false) to disable recording.

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

If no any other operation, you can use pass `nullptr` and it will automatically be replaced by an empty lambda.
   
**All these function pointers can be lambda. When constructed, all these function pointers has an empty lambda as default value.**

### `GeneIt_t`
It's same as iterator type of stl list containing `Gene` types.

<br>

## Member details
### `_eliteIt`
Iterator to elite gene.

### `_population`
A list to store all Genes.

### `_option`
Options of solver.

### `_generation`
Generations used.

### `_failTimes`
Generations that failed to updated elite.

### `_record`
Vector to store elite fitness history in each generation.

**Exists only when template parameter [`Record`](#record) is true.**

### `_args`
Other parameters pack.

### `_initializeFun`
Function pointer to initialize.

### `_fitnessFun`
Function pointer to calculate fitness value.

### `_crossoverFun`
Function pointer to apply crossover between 2 `Var_t` and produce 2 more `Var_t`.

### `_mutateFun`
Function pointer to apply mutation of a `Var_t`.

### `_otherOptFun`
Function pointer to apply special other operation after each generation.

<br>

## Function details
### `GABase()`
Default constructor. All function pointers are initialized with empty lambda.

### `~GABase()`
Default virtual constructor.

### `initialize()`
Initiaze solver with all 5 function pointers, other parameters and `GAOption`.

### `run()`
Start solving. Must be called after initialization.

### `eliteIt() const`
Get iterator to result gene.

### `result() const`
Get result `Var_t`.

### `population() const`
Get the whole population.

### `record() const`
Get record vector.

**Exists only when template parameter [`Record`](#record) is true.**

### `option() const`
Get `GAOption`.

### `generation() const`
Get generation used.

### `failTimes() const`
Get generations that failed to updated elite.

### `args() const`
Get parameter pack.

### `calculateAll()`
Calculate fitness values for every genes.

This is a virtual function, so you can cover it by inheriting.

### `select()=0`
Apply selection to select better genes.

This is a pure virtual function since comparsion between abstrcat `Fitness_t` is impossible -- you don't even know how many object functions there are!

### `crossover()`
Apply crossover operation which let randomly 2 gene born 2 more new one. They will be added to population.

### `mutate()`
Apply mutation operation which slightly modify gene.
