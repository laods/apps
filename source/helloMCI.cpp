#include <iostream>

#include <opm/core/utility/MonotCubicInterpolator.hpp>

using namespace Opm;
using namespace std;

int main(int varnum, char** vararg)
{
  vector<double> x(3, 0.0);
  vector<double> f(x);
  x[1] = 1.0;
  x[2] = 2.0;
  f[1] = 1.0;
  f[2] = 4.0;
  MonotCubicInterpolator mci(x, f);
  cout << "Hello world\n";

  return 0;
    
}
