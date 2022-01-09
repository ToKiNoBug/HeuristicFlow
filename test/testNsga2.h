#ifndef TESTNSGA2_H
#define TESTNSGA2_H

#ifdef OptimT_NO_OUTPUT
#undef OptimT_NO_OUTPUT
#endif

#include <OptimTemplates/Genetic>
#include <array>

//test with 2 object function
//object1: f1=4x^2+4y^2
//object2: f2=(x-5)^2+(y-5)^2
//0<=x<=5,0<=y<=3
class testNsga2 :
        public OptimT::GABase<
        std::array<double,2>,//Var_t
        std::array<double,2>,//Fitness_t
        OptimT::RECORD_FITNESS>
{
public:
    testNsga2();
    using Base_t =
        GABase<std::array<double,2>,std::array<double,2>,OptimT::RECORD_FITNESS>;
    class GeneItPlus_t
    {
    public:
        GeneItPlus_t() {};
        ~GeneItPlus_t() {};
        Base_t::GeneIt_t iterator;
        uint32_t idx;
        int paretoRank;
        bool isSelected;
        std::array<double,2> congestion;
    };

    void paretoFront(std::vector<const Base_t::Gene*>&) const;

protected:
    static bool isBetter(const std::array<double,2>&,const std::array<double,2>&);
    void select();
};

void runNSGA2();

#endif // TESTNSGA2_H
