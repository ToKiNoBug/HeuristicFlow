/*
 Copyright © 2021-2022  TokiNoBug
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

#include <Eigen/Dense>
#include <HeuristicFlow/Genetic>
#include <iostream>
#include <ctime>
using namespace Eigen;
using namespace std;

// Test SOGA with TSP problem
void testTSP_SOGA(const uint32_t PointNum) {
  static const uint8_t DIM = 2;
  // static const double LengthBase=100;
  typedef array<double, DIM> Point_t;
  vector<Point_t> points;
  points.clear();
  points.reserve(PointNum);
  for (uint32_t i = 0; i < PointNum; i++) {
    points.emplace_back();
    for (auto &v : points.back()) {
      v = heu::randD();
    }
  }
  /*
  cout<<"Generated random points :"<<endl;
  for(const auto & i : points) {
      cout<<'(';
      for(auto j : i) {
          cout<<j<<" , ";
      }
      cout<<")"<<endl;
  }
  */
  typedef pair<double, uint32_t> permUnit;

  // typedef vector<permUnit> permulation;
  //       var,        less=better,    data src
  heu::SOGA<vector<double>, heu::FITNESS_LESS_BETTER, heu::DONT_RECORD_FITNESS, heu::Tournament,
            std::tuple<const vector<Point_t> *>>
      algo;
  static const uint8_t dataIdx = 0;
  typedef tuple<const vector<Point_t> *> Args_t;
  Args_t args;  //=make_tuple(PointNum,points.data());
  get<dataIdx>(args) = &points;
  // initialize function

  auto initializeFun = [](vector<double> *x, const Args_t *_args) {
    const auto permL = get<dataIdx>(*_args)->size();
    x->resize(permL);
    for (double &var : *x) var = heu::randD();
  };

  // calculate fitness
  auto calculateFun = [](const vector<double> *x, const Args_t *_args, double *f) {
    const auto permL = x->size();
    vector<permUnit> perm(permL);
    for (uint32_t i = 0; i < permL; i++) {
      perm[i].first = x->operator[](i);
      perm[i].second = i;
    }
    std::sort(perm.begin(), perm.end(),
              // compare function for std::sort
              [](const permUnit &a, const permUnit &b) { return a.first > b.first; });

    double L = 0;
    for (uint32_t i = 1; i < permL; i++) {
      const Point_t &prev = get<dataIdx>(*_args)->operator[](perm[i - 1].second);
      const Point_t &cur = get<dataIdx>(*_args)->operator[](perm[i].second);
      double curL = 0;
      for (uint8_t d = 0; d < DIM; d++) {
        curL += (prev[d] - cur[d]) * (prev[d] - cur[d]);
      }
      L += curL;
    }
    *f = L;
  };

  // calculate natural perm pathL
  {
    vector<double> natural(PointNum);
    for (uint32_t i = 0; i < PointNum; i++) {
      natural[i] = double(i) / PointNum;
    }
    double naturalPathL;
    calculateFun(&natural, &args, &naturalPathL);
    cout << "default pathL = " << naturalPathL << endl;
  }

  auto crossoverFun = heu::GADefaults<vector<double>, Args_t>::cFunXd<heu::DivCode::DivCode_Half>;

  auto mutateFun = [](const vector<double> *src, vector<double> *x, const Args_t *) {
    *x = *src;
    const auto permL = x->size();
    if (heu::randD() < 0.5) {  // modify one element's value
      double &v = x->operator[](std::rand() % permL);
      v += heu::randD(-0.01, 0.01);
      v = std::min(v, 1.0);
      v = std::max(v, 0.0);
    } else {
      const uint32_t flipB = std::rand() % (permL - 1);
      const uint32_t flipE = std::rand() % (permL - flipB - 1) + flipB + 1;
      // cout<<"flipB="<<flipB<<" , flipE"<<flipE<<endl;
      const uint32_t flipN = flipE - flipB + 1;
      for (uint32_t flipped = 0; flipped * 2 <= flipN; flipped++) {
        swap(x->operator[](flipB + flipped), x->operator[](flipE - flipped));
      }
    }
  };

  heu::GAOption opt;
  opt.maxGenerations = 30 * PointNum;
  opt.maxFailTimes = opt.maxGenerations;
  //  If maxFailTimes is equal or greater than maxGenerations, the solver will definatly runs
  //  for 100 generations. Since member `maxFailTimes` is of type `size_t`, you can also use
  //  -1 to achieve this, however it's kind of undefined behavior

  algo.setiFun(initializeFun);
  algo.setmFun(mutateFun);
  algo.setfFun(calculateFun);
  algo.setcFun(crossoverFun);
  algo.setOption(opt);
  algo.setArgs(args);
  algo.initializePop();

  cout << "run!" << endl;
  std::clock_t c = std::clock();
  algo.run();
  c = std::clock() - c;
  cout << "finished with " << algo.generation() << " generations and " << double(c) / CLOCKS_PER_SEC
       << " s\n";
  cout << "result fitness = " << algo.bestFitness() << endl;
}

int main(int argc,char**argv) {
  bool is_auto=false;

  for(int i=0;i<argc;i++) {
    if(std::string_view{argv[i]}=="--auto") {
      is_auto=true;
      break;
    }
  }

  int NodeNum = 100;
  if(!is_auto) {
    cout << "Input node number : ";
    cin >> NodeNum;
  }
  testTSP_SOGA(NodeNum);
  // system("pause");
  return 0;
}
