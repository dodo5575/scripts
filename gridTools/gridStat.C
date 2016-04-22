// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"

using namespace std;

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 2 ) {
    printf("Usage: %s srcGrid\n", argv[0]);
    return 0;
  }

  const char* srcGrid = argv[1];  
  Grid src(srcGrid);
  
  // Get extreme values.
  int minInd = 0;
  double minVal = src.getValue(0);
  int maxInd = 0;
  double maxVal = src.getValue(0);

  // const double level = 0.0;
  const int n = src.length();
  double sum = 0.0;
  double sumSq = 0.0;

  Vector3 lowSum(0.0);
  double lowWeight = 0.0;

  for (int i = 0; i < n; i++) {
    double v = src.getValue(i);
    Vector3 r = src.getPosition(i);
    sum += v;
    sumSq += v*v;

    lowSum += exp(-v)*r;
    lowWeight += exp(-v);

    if (v < minVal) {
      minVal = v;
      minInd = i;
    }
    
    if (v > maxVal) {
      maxVal = v;
      maxInd = i;
    }
  }

  double mean = sum/n;
  double std = sqrt((sumSq - sum*sum/n)/(n-1));
  printf("%g %g\n", sumSq, sum*sum);
  
  printf("count: %d\n", n);
  printf("mean: %.12g\n", mean);
  printf("stdev: %.12g\n", std);
  
  Vector3 minPos = src.getPosition(minInd);
  Vector3 maxPos = src.getPosition(maxInd);
  double boltzMean = lowWeight/n;
  Vector3 boltzCen = lowSum/lowWeight;
  printf("min: %.12g %.12g %.12g %.12g\n", minVal, minPos.x, minPos.y, minPos.z);
  printf("max: %.12g %.12g %.12g %.12g\n", maxVal, maxPos.x, maxPos.y, maxPos.z);
  printf("boltzmann: %.12g %.12g %.12g %.12g\n", boltzMean, boltzCen.x, boltzCen.y, boltzCen.z);

  return 0;
}
