
#include <Eigen/Dense>
#include <HeuristicFlow/PSO>
#include <cmath>
#include <iostream>
using namespace Eigen;
using namespace std;

// this tests shows the hierarchy of PSO solvers.
void testPSOBase() {
  using Var_t = std::vector<double>;

  heu::internal::PSOAbstract<Var_t, double, heu::DONT_RECORD_FITNESS, void, nullptr, nullptr>* abstractNoRec = nullptr;

  heu::internal::PSOAbstract<Var_t, double, heu::RecordOption::RECORD_FITNESS, void, nullptr, nullptr>* abstractDoRec =
      nullptr;

  heu::internal::PSOBase<Var_t, 10, double, heu::RecordOption::DONT_RECORD_FITNESS, void, nullptr, nullptr>* baseNoRec =
      nullptr;

  heu::internal::PSOBase<Var_t, 0, double, heu::RecordOption::RECORD_FITNESS, void, nullptr, nullptr>* baseDoRec =
      nullptr;

  // base is derived from abstract
  abstractNoRec = baseNoRec;
  abstractDoRec = baseDoRec;

  // DoRec is derived from NoRec
  abstractNoRec = abstractDoRec;
}

int main() {
  testPSOBase();
  system("pause");
  return 0;
}
