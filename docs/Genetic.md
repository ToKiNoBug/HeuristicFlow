# Genetic algorithm template class

GA is an templated implementation of genetic algorithm using elitism strategy. 

<br>
<br>

## Defination
```cpp
namespace OptimT
{
template<typename Var,bool isGreaterBetter,class ...Args> 
class GA;
}
```

In the template above, `Var` is the variable to be modified, while `isGreaterBetter` denotes that whether a greater fitness value means a better variable. `...Args` means extra custom parameters, they can be some datas to calculate the fitness of a variable like trainning data set, or some other parameters like maximum and minimum value of `Var`, or even trainning rate. Also you can leave it blank. Several examples:
```cpp
OptimT::GA<std::array<double,16>,true> solver_empty;
OptimT::GA<BPNetwork_t,false,DataSet_t> GABP_trainner;
OptimT::GA<Eigen::Array4d,true,std::vector<Point_t>,double,double,double> TSP_solver;
```

<br>
<br>

## General Options
In namespace `OptimT` I defined a non-template struct `GAOptions`, including several hyper-parameters of all genetic algorithm will use. They are listed as below:

| Variable name | Type | Defaule value | Description |
| :----: | :----: | :----: | :---- |
| populationSize | `size_t` | 100 | Size of population. |
| maxGeneration | `size_t` | 300 | Maximum generation. GA will stop once reached this limitation. |
| maxFailTimes | `size_t` | 50 | GA will stop once best solve hasn't been improved for continuous `maxFailTimes` generations. |
| crossoverProb | `double` | 0.8 | Probability of a non-elite individual to join crossover. |
| mutateProb | `double` | 0.05 | Probability of a non-elite individual to get mutated. |

All these members are public in order to be modified conveniently.

<br>
<br>

## Types inside GA
1. `ArgsType`
    ```cpp
    typedef std::tuple<Args...> ArgsType;
    ```
    In class `GA`, `...Args` will be packed in a tuple (`std::tuple`) and such tuple are called **args**, which means parameters of this optimization task.
2. `initializeFun`
    ```cpp
    typedef void(* initializeFun)(Var*,const ArgsType*);
    ```
    Function pointer type to initialize an individual in population. This function will be called only when initializing. For best flexibility, a constant pointer to parameters is provided.
    
3. `fitnessFun`
   ```cpp
    typedef double(* fitnessFun)(const Var*,const ArgsType*);
   ```
   Fucntion pointer type to calculate fitness for any individual. Although it may be not too flexible, I assign that fitness value are `double`.

   This function will be called to each individual in each generation, so it often spend most time when optimizing. So when you are solving a complex problem, make this function as fast as possible.

4. `crossoverFun`
   ```cpp
   typedef void(* crossoverFun)(const Var*,const Var*,Var*,Var*,const ArgsType*);
   ```
   Function pointer type to do crossover. The first 2 constant pointers can be called *parent1* and *parent2* while call rest non-constant pointers *child1* and *child2*. When doing crossover, parents are chosen from population while children are newly inserted to population. This function should assign children's value by parents' value.
   Note that children are constructed with **default constructor**, **make sure that your `Var` type does have a default constructor**.

5. `mutateFun`
    ```cpp
    typedef initializeFun mutateFun;
    ```
    Function pointer type to do mutation. It has same type with `initializeFun` but they have different meanings. This function should slightly modify given `Var`. Same as above, these functions have constant access to parameters.

6. class `Gene`
   ```cpp
   GA<Var,bool,...>::class Gene
   ```
   Gene type encapsulating a `Var` object. It has following members (and functions):
   | Type/Return Type | Name | Type |
   | ----: | :---- | :----: |
   | `Var` | `self` | public member |
   | `void` | `setUncalculated()` | public function |
   | `bool`  | `isCalculated() const` | public function |
   | `double` | `fitness() const` | public function |
   | `bool` | `_isCalculated` | protected member |
   | `double` | `_Fitness` | protected member |

    To avoid repeated fitness caculation, an extra boolean is added to note whether a `Gene` has already get its fitness value. You can call `void Gene::setUncalculated()` to set a Gene to uncalculated, however, only class `GA` is able to assign fitness to `Gene` and turn it into calculated. This can be done since `GA` is a friend class of `Gene`.
7. `otherOptFun`
   ```cpp
   typedef void(* otherOptFun)
        (ArgsType*,std::list<Gene>*,size_t generation,size_t failTimes,const GAOption*);
   ```
   It means "other operate function". This function is called in every generation, after selection but before crossover. In this function, you can modify not only args but also the whole population. In most condition this step is not needed, but I love flexibility.
   
   **All these function pointers can be lambda. When constructed, all these function pointers has an empty lambda as default value.**

8. `GeneIt`(**protected**)
   ```cpp
    typedef typename std::list<Gene>::iterator GeneIt;
    ```
    It's same as iterator type of stl list containing `Gene` types. GA uses std::list to store the whole population since it's fast for inserting and erasing while random access is merely needed.

## Public functions of GA
| Type | Defination | Description |
| ----: | :---- | :---- |
|  | `GA()` | Constructor |
| `virtual` | `~GA()` | Destructor |
| `virtual void` | `initialize(initializeFun,fitnessFun,crossoverFun,mutateFun,otherOptFun=nullptr,const GAOption&=GAOption(),const ArgsTyp&=ArgsType())` | initializer |
| `virtual void` | `run()` | Run algorithm |
| `typename std::list<Gene>::iterator` | `eliteIt() const` | Get iterator to elite |
| `const Var &` | `result() const` | Get result Var |
| `const std::list<Gene> & ` | `population() const` | Get whole population |
| `const std::vector<double> & ` | `recording() const` | Get trainning curve |
| `const GAOption & ` | `option() const` | Get `GAOption` member |
| `size_t `| `generation() const` | Get generations used |
| `size_t ` | `failTimes() const` | Get fail times used |
| `const ArgsType & ` | `args() const` | Get parameter pack(`std::tuple`) |
 
GA don't have any public member.