/*
 Copyright © 2022  TokiNoBug
This file is part of OptimTemplates.

    OptimTemplates is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OptimTemplates is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OptimTemplates.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef ENUMERATIONS_HPP
#define ENUMERATIONS_HPP

#include <stdint.h>

namespace OptimT {
    
///whether to record trainning curve of not
enum RecordOption : uint8_t {
    RECORD_FITNESS=true,
    DONT_RECORD_FITNESS=false
};

///convert enumeration to string
inline const char * Enum2String(RecordOption r) {
    switch (r) 
    {
        case RECORD_FITNESS:
        return "RECORD_FITNESS";
        case DONT_RECORD_FITNESS:
        return "DONT_RECORD_FITNESS";
    }
}

///optimization direction
enum FitnessOption : uint8_t {
    FITNESS_LESS_BETTER=false,
    FITNESS_GREATER_BETTER=true,
};

///convert enumeration to string
inline const char * Enum2String(FitnessOption f) {
    switch (f)
    {
        case FITNESS_LESS_BETTER:
        return "FITNESS_LESS_BETTER";
        case FITNESS_GREATER_BETTER:
        return "FITNESS_GREATER_BETTER";
    }
}

///whether it's constrainted
enum ConstraintOption : uint8_t {
    NONCONSTRAINT,
    IS_CONSTRAINT
};
///convert enumeration to string
inline const char * Enum2String(ConstraintOption c) {
    switch (c)
    {
        case NONCONSTRAINT:
        return "NONCONSTRAINT";
        case IS_CONSTRAINT:
        return "IS_CONSTRAINT";
    }
}

///which type of vector to use
enum DoubleVectorOption {
    Std='S',
    Eigen='E',
    Custom='C'
};
///convert enumeration to string
inline const char * Enum2String(DoubleVectorOption e) {
    switch (e) 
    {
        case Std:
        return "C++ std vector/array";
        case Eigen:
        return "Eigen Array";
        default:
        return "Unknown custom types";
    }
}

}

#endif // ENUMERATIONS_HPP
