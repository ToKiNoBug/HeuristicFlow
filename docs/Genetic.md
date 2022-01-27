# Genetic algorithm template class

Genetic has some implementations of genetic algorithm solver.

<br>
<br>

## Structs/Classes
1. [GAOption](./Genetic/GAOption.md)
2. [GABase](./Genetic/GABase.md)
3. [SOGA](./Genetic/SOGA.md)
4. MOGAAbstract
5. MOGABase
6. NSGA2Base
7. [NSGA2](./Genetic/NSGA2.md)

## Class diagram
```mermaid
classDiagram

class GAOption 

class GABase {
    <<abstract>>

}

class SOGA 

class MOGAAbstract {
    <<abstract>>
}

class MOGABase {
    <<abstract>>
}

class NSGA2Base {
    <<abstract>>
}

class NSGA2

GAOption *-- GABase
GABase <|-- SOGA
GABase <|-- MOGAAbstract
MOGAAbstract <|-- MOGABase
MOGABase <|-- NSGA2Base
NSGA2Base <|-- NSGA2

```