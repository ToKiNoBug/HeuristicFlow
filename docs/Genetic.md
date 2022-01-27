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

class GAOption {
    public size_t populationSize
    public size_t maxGeneration
    public size_t maxFailTimes
    public double crossoverProb
    public double mutationProb
}

class GABase {
    <<abstract>>
    protected GAOption _opt
    protected iFun_t iFun
    protected fFun_t fFun
    protected cFun_t cFun
    protected mFun_t mFun
    protected ooFun_t ooFun
    public virtual void initialize();
    public virtual void run();
    protected virtual void select()*;
    protected virtual void crossover();
    protected virtual void mutate()*;

}

class SOGA {
    protected GeneIt_t _elite
    protected virtual void select();
    protected virtual void mutate();
}

class MOGAAbstract {
    <<abstract>>
    public void paretoFront()
    public const auto & pfGenes()const;
    protected unordered_map _pfGenes
}

class MOGABase {
    <<abstract>>
}

class NSGA2Base {
    <<abstract>>
}

GAOption o-- GABase
GABase <|-- SOGA
GABase <|-- MOGAAbstract
MOGAAbstract <|-- MOGABase
MOGABase <|-- NSGA2Base
NSGA2Base <|-- NSGA2

```