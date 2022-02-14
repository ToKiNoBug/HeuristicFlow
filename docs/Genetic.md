# Genetic algorithm template class

This module has some implementations of varities of genetic algorithm solvers, including `SOGA` for single objectives while `NSGA2` and `NSGA3` for multiple objectives.

<br>
<br>

## Structs/Classes
1. [GAOption](./Genetic/GAOption.md)
2. [GABase](./Genetic/GABase.md)
3. [SOGA](./Genetic/SOGA.md)
4. [MOGAAbstract](./Genetic/MOGAAbstract.md)
5. [MOGABase](./Genetic/MOGABase.md)
6. [NSGA2Base](./Genetic/NSGA2Base.md)
7. [NSGA2](./Genetic/NSGA2.md)

## Class diagram

Class name in italic refers to an abstract class, while that in bold means that the class is designed to be used directly. 2 classes in a dashline box means a parital specialization for different conditions, such as to record traninning curve or not.

![class diagram in image](https://raw.githubusercontent.com/ToKiNoBug/SlopeCraftTutorial/Images4HeuristicFlow/Genetic/classDiagram.png)

Since there're so many unrelated options in template parameters(record fitness or not, dynamic objectives or not, std vectors or eigen vectors), such a complex structure is applied to seperate each choice and to improve code reusability. 

In each inheritance, the template matching principle ensures that all specializations will be used correctly, while the polymorphism mechanic enables us to replace previously defined virtual functions perfectly which are called in base class with our new one. 

It's really complicated but at least such strcuture freed me from having to copy same code for multiple times and apply every single updation to every copies of them.