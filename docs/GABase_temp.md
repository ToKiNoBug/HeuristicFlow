# Genetic algorithm template class

GA is an templated implementation of genetic algorithm using elitism strategy. 

<br>
<br>

## Defination
```cpp
namespace OptimT
{
template<typename Var_t,bool isGreaterBetter,class ...Args> 
class GA;
}
```

In the template above, `Var_t` is the Var_tiable to be modified, while `isGreaterBetter` denotes that whether a greater fitness value means a better Var_tiable. `...Args` 

<br>
<br>

## General Options
In namespace `OptimT` I defined a non-template struct `GAOptions`, including several hyper-parameters of all genetic algorithm will use. They are listed as below:

| Var_tiable name | Type | Defaule value | Description |
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
2. `initializeFun`
    ```cpp
    typedef void(* initializeFun)(Var_t*,const ArgsType*);
    ```
    
3. `fitnessFun`
   ```cpp
    typedef double(* fitnessFun)(const Var_t*,const ArgsType*);
   ```

4. `crossoverFun`
   ```cpp
   typedef void(* crossoverFun)(const Var_t*,const Var_t*,Var_t*,Var_t*,const ArgsType*);
   ```

5. `mutateFun`
    ```cpp
    typedef initializeFun mutateFun;
    ```
    

6. class `Gene`
   ```cpp
   GA<Var_t,bool,...>::class Gene
   ```
   
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
| `const Var_t &` | `result() const` | Get result Var_t |
| `const std::list<Gene> & ` | `population() const` | Get whole population |
| `const std::vector<double> & ` | `recording() const` | Get trainning curve |
| `const GAOption & ` | `option() const` | Get `GAOption` member |
| `size_t `| `generation() const` | Get generations used |
| `size_t ` | `failTimes() const` | Get fail times used |
| `const ArgsType & ` | `args() const` | Get parameter pack(`std::tuple`) |
 
GA don't have any public member.