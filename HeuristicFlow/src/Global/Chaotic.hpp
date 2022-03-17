// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifndef Heu_CHAOTIC_HPP
#define Heu_CHAOTIC_HPP
#include "./LogisticChaos.hpp"
namespace Heu {

extern LogisticChaos global_logistic;

///logistic random number in range (0,1)
inline double logisticD() {
    return global_logistic();
}

}   // Heu

#endif // CHAOTIC_HPP