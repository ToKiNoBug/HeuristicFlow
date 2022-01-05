# OptimT::GAOption
Defined in header [GABase.h](../../GA/GABase.h)

<br>

## Defination
```cpp
    struct GAOption
```
`GAOption` is a struct encapsulating several general options for most types of genetic algorithm. It's a simple struct with all members public and without any function except constructor.

<br>

## Member types
| Type | Name | Default value | Description |
| ----: | :---- | :----: | :---- |
| `size_t` | `populationSize` | 100 | Size of population. |
| `size_t` | `maxGeneration` | 300 | Maximum generation. GA will stop once reached this limitation. |
| `size_t` | `maxFailTimes` | 50 | GA will stop once best solve hasn't been improved for continuous `maxFailTimes` generations. |
| `double` | `crossoverProb` | 0.8 | Probability of a non-elite individual to join crossover. |
| `double` | `mutateProb` | 0.05 | Probability of a non-elite individual to get mutated. |


## Member functions
| Return type | Defination | Description |
| ----: | :---- | :---- |
| (Constructor) | `GAOption()` | Construct and initialize all members to default value |