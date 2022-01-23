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

#include <iostream>
#include <ctime>
#include "def_TestingPSO.h"

#include <unordered_set>

using namespace std;

OptimT_MAKE_GLOBAL

int main()
{
    testTSP(500);
    /*
    static const size_t N=10000000;
    static const size_t POW=10;
    Eigen::ArrayXd a,b;
    a.setLinSpaced(N,1,2);
    system("pause");
    clock_t c=clock();
    for(size_t i=0;i<N;i++) {
        a(i)=
        OptimT::power<POW>(a(i))
        //OptimT::power(a(i),POW)
        ;
    }
    c=clock()-c;
    cout<<"Mean time cost = "<<1e9*double(c)/CLOCKS_PER_SEC/N<<" nano seconds"<<endl;
    cout<<a.mean()<<endl;
    */
    //cout<<OptimT::power(1.01,890)<<endl;
    system("pause");
    return 0;
}
