// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef EIGEN_HEU_LOGISTICCHAOS_HPP
#define EIGEN_HEU_LOGISTICCHAOS_HPP

#include <cstdint>

namespace Eigen
{

class LogisticChaos
{
  public:
    LogisticChaos()
    {
      _value=0.1;
    }

    LogisticChaos(double seed)
    {
      _value=seed;
    }

    double operator()()
    {
      _value*=miu*(1-_value);
      return _value;
    }

    double value() const
    {
      return _value;
    }

    void makeSequence(double *dst,size_t L)
    {
      if(L<=0) return;
        dst[0]=operator()();
      for(size_t i=1;i<L;i++)
      {
        dst[i]=miu*dst[i-1]*(1-dst[i-1]);
      }
      _value=dst[L-1];
    }

    void iterate(size_t It)
    {
      for(size_t i=0;i<It;i++)
      {
        _value*=miu*(1-_value);
      }
    }
    
  private:
    double _value;
    constexpr static const double miu=4.0;
};

}
#endif // EIGEN_HEU_LOGISTICCHAOS_H
