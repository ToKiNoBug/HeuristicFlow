// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifndef Heu_CHAOSBASE_HPP
#define Heu_CHAOSBASE_HPP

#include <type_traits>

namespace Heu {

template<typename ele_ts>
class ChaoseBase
{
protected:
    ChaoseBase() {};
    ~ChaoseBase() {};
};

};  //  namespace Heu

#endif  //  Heu_CHAOSBASE_HPP
