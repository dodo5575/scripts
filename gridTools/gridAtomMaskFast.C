// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "BaseGrid.H"
#include "Grid.H"
#include "CellDecomposition.H"
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

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 6 ) {
    printf("Make a mask set to 1 for grid nodes within radius of any point.\n");
    printf("Usage: %s srcGrid pointFile radius blurCount outFile\n", argv[0]);
    return 0;
  }

  const char* srcFile = argv[1];
  const char* pointFile = argv[2];
  const double radius = strtod(argv[3], NULL);
  const int blurCount = atoi(argv[4]);
  const char* outFile = argv[argc-1];

  // Load the grids.
  Grid src(srcFile);
  src.zero();
  int n = src.length();
  printf("Loaded `%s' which contains %d nodes.\n", srcFile, n); 

  // Load the coordinates.
  printf("Loading the coordinates.\n");
  Scatter pos(pointFile);
  printf("Loaded %d points from `%s'.\n", pos.length(), pointFile);

  // Generate the cell decomposition.
  double cutoff = 1.2*radius;
  Matrix3 sys = src.getBox();
  Vector3 origin = src.getOrigin();
  printf("\nGenerating cell decomposition for %s with cutoff = %.10g...\n", argv[1], cutoff);
  CellDecomposition cell(sys, origin, cutoff);
  printf("Cell decomposition of %i cells\n", cell.length());
  cell.decompose(pos);
  printf("Cell decomposition contains %i points\n", cell.countPoints());

  // Make the mask.
  Grid mask(src);
  int complete = 0;
  for (int i = 0; i < n; i++) {
    Vector3 r = src.getPosition(i);
    // Get the indices (in pos) of possible neighboring points.
    IndexList neigh = cell.neighborhood(r);
    const int nNeighs = neigh.length();

    for (int j = 0; j < nNeighs; j++) {
      // Get the position of each possible neighboring point.
      Vector3 p = pos.get(neigh.get(j));
      if (within(p, r, radius, src)) {
	mask.setValue(i, 1.0);
	break;
      }
    }

    // Write the progress.
    int comp = (100*i)/n;
    if (abs(complete - comp) >= 5) {
      printf("%d percent complete\n", comp);
      complete = comp;
    }
  }

  // Blur.
  for (int j = 0; j < blurCount; j++) mask.blur();

  char comments[256];
  sprintf(comments, "%s patched", srcFile);
  mask.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
