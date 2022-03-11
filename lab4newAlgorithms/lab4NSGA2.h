/*
 Copyright © 2022  TokiNoBug
This file is part of HeuristicFlow.

    HeuristicFlow is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HeuristicFlow is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HeuristicFlow.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef TESTNSGA2_H
#define TESTNSGA2_H

#include "includes.h"
#include <array>
#include <iostream>
//test with 2 object function
//object1: f1=4x^2+4y^2
//object2: f2=(x-5)^2+(y-5)^2
//0<=x<=5,0<=y<=3
class testNsga2 :
        public Heu::GABase<
        std::array<double,2>,//Var_t
        std::array<double,2>,//Fitness_t
        Heu::RECORD_FITNESS,void,nullptr,nullptr,nullptr,nullptr>
{
public:
    testNsga2();
    using Base_t =
        Heu::GABase<std::array<double,2>,std::array<double,2>,Heu::RECORD_FITNESS,
            void,nullptr,nullptr,nullptr,nullptr>;
    struct infoUnit
    {
    public:
        //bool isLayered;
        bool isSelected;
        int32_t sortType;
        uint32_t index;
        uint32_t domainedByNum;
        Base_t::GeneIt_t iterator;
        std::array<double,2> congestion;
    };

    void paretoFront(std::vector<const Base_t::Gene*>&) const;

    std::array<double,2> bestFitness() const;

protected:
    static bool isBetter(const std::array<double,2>&,const std::array<double,2>&);
    void select();

    void mutate();
};

void runNSGA2();

#endif // TESTNSGA2_H
