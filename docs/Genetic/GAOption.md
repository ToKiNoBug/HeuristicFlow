# OptimT::GAOption
`GAOption` is a struct encapsulating several general options for most types of genetic algorithm.

| Header: | `#include<Genetic>` |
| ----: | :---- |
| Location: | [GAOption.hpp](../../GA/GAOption.hpp) |

<br>

## Defination
```cpp
    struct OptimT::GAOption;
```

<br>

## Members
| Access | Type | Name | Default value |
| :----: | ----: | :---- | :----: |
| public | `size_t` | [`populationSize`](#populationsize) | 100 |
| public | `size_t` | [`maxGeneration`](#maxgeneration) | 300 |
| public | `size_t` | [`maxFailTimes`](#maxfailtimes) | 50 |
| public | `double` | [`crossoverProb`](#crossoverprob) | 0.8 |
| public | `double` | [`mutateProb`](#mutateprob) | 0.05 |

<br>

## Member functions
| Access | Return type | Defination |
| :----: | ----: | :---- |
| public |  | [`GAOption()`](#gaoption) |

<br>

## Detailed description
`GAOption` is a non-template struct, which means that it's suitable to most genetic algorithm. It's a simple struct with all members public and without any function except constructor.

<br>

## Member details
### `populationSize`
Size of population. 

### `maxGeneration`
Maximum generation. GA will stop once reached this limitation.

### `maxFailTimes`
GA will stop once best solve hasn't been improved for continuous `maxFailTimes` generations.

### `crossoverProb`
Probability of a non-elite individual to join crossover.

### `mutateProb`
Probability of a non-elite individual to get mutated.

<br>

## Function details
### `GAOption()`
Construct and initialize all members to their default values.
