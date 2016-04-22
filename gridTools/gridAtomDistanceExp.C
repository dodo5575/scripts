// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"
#include "Scatter.H"

using namespace std;

bool within(const Scatter& pos, Vector3 r, double rMax) {
  const int n = pos.length();
  double r2 = rMax*rMax;
  for (int p = 0; p < n; p++) {
    Vector3 d = pos.get(p) - r;
    if (d.length2() < r2) return true;
  }
  return false;
}

bool within(const Scatter& pos, Vector3 r, double rMax, const Grid& g) {
  const int n = pos.length();
  double r2 = rMax*rMax;
  for (int p = 0; p < n; p++) {
    Vector3 d = g.wrapDiff(pos.get(p) - r);
    if (d.length2() < r2) return true;
  }
  return false;
}

bool within(Vector3 p, Vector3 r, double rMax, const Grid& g) {
  double r2 = rMax*rMax;
  Vector3 d = g.wrapDiff(p - r);
  return (d.length2() < r2);
}

// Calculate a + b*x^c.
double powerLaw(double x, double a, double b, double c) {
  return a + b*pow(x,c);
}

// Calculate a + b*exp(c*x).
double exponential(double x, double a, double b, double c) {
  return a + b*exp(c*x);
}


///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 8 ) {
    printf("Set the grid v = a + b*exp(c*D), where d is the D distance to the nearest atom.\n");
    printf("Usage: %s srcGrid pointFile a b c outFile\n", argv[0]);
    return 0;
  }

  const char* srcFile = argv[1];
  const char* pointFile = argv[2];
  const double a = strtod(argv[3], NULL);
  const double b = strtod(argv[4], NULL);
  const double c = strtod(argv[5], NULL);
  const double cutDist = strtod(argv[6], NULL);
  const char* outFile = argv[argc-1];

  // Load the grids.
  Grid src(srcFile);
  src.zero();
  int n = src.length();
  printf("Loaded `%s' which contains %d nodes.\n", srcFile, n); 

  // Load the coordinates.
  printf("Loading the coordinates.\n");
  Scatter pos(pointFile);
  int posNum = pos.length();
  printf("Loaded %d points from `%s'.\n", pos.length(), pointFile);

  int complete = 0;
  for (int i = 0; i < n; i++) {
    Vector3 r = src.getPosition(i);

    Vector3 dr = pos.get(0)-r;
    double minDist = dr.length();
    for (int j = 1; j < posNum; j++) {
      Vector3 dr = pos.get(j)-r;
      double dist = dr.length();

      if (dist < minDist) minDist = dist;
    }
    if (minDist < cutDist) minDist = cutDist;

    double v = exponential(minDist, a, b, c);
    src.setValue(i, v);

    // Write the progress.
    int comp = (100*i)/n;
    if (abs(complete - comp) >= 5) {
      printf("%d percent complete\n", comp);
      complete = comp;
    }
  }
  char comments[256];
  sprintf(comments, "%s distance map", srcFile);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
