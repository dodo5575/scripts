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

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 7 ) {
    printf("Match one grid to another in the region between radiusNear and radiusFar\n");
    printf("Paste the patch up to radiusFar.\n");
    printf("Usage: %s srcGrid patchGrid coordinateFile radiusNear radiusFar outFile\n", argv[0]);
    return 0;
  }

  const char* srcFile = argv[1];
  const char* patchFile = argv[2];
  const char* coorFile = argv[3];
  const double radNear = strtod(argv[4], NULL);
  const double radFar = strtod(argv[5], NULL);
  const char* outFile = argv[argc-1];


  // Load the grids.
  Grid src(srcFile);
  printf("Loaded `%s'.\n", srcFile); 
  int n = src.length();
  Grid patch(patchFile);
  printf("Loaded `%s'.\n", patchFile); 
  
  // Determine the number of blurs.
  double blurCount = 6;

  // Load the coordinates.
  printf("Loading the coordinates. ");
  Scatter pos(coorFile);

  printf("Loaded %d points from `%s'.\n", pos.length(), coorFile);

  // Make masks based on the coordinate file.
  printf("Generating the masks.\n");
  Grid maskFar(src);
  maskFar.zero();
  Grid maskNear(maskFar);
  for (int i = 0; i < n; i++) {
    Vector3 r = src.getPosition(i);
    if (within(pos, r, radNear)) maskNear.setValue(i, 1.0);
    if (within(pos, r, radFar)) maskFar.setValue(i, 1.0);
  }

  // Smooth the masks.
  printf("Blurring the masks.\n");
  for (int b = 0; b < blurCount; b++) {
    maskFar.blur();
    maskNear.blur();
  }

  // Make a mask for the region between far and near.
  Grid maskBorder(maskFar);
  maskBorder.alphaDifference(maskNear);

  // Write the masks.
  //maskFar.write("maskFar.dx");
  //maskNear.write("maskNear.dx");
  //maskBorder.write("maskBorder.dx");
  char maskName[256];
  sprintf(maskName, "%s.mask", outFile);
  maskFar.write(maskName, "far mask");
  sprintf(maskName, "%s.mask", outFile);
  maskBorder.write(maskName, "border mask");
  
  // Get the average.
  double srcV = src.averageRegion(maskBorder);
  double patchV = patch.averageRegion(maskBorder);
  printf("Average values: %g %g\n", srcV, patchV);
  
  // Shift the potential.
  patch.shift(srcV-patchV);   

  // Paste the patch.
  src.paste(patch, maskFar);

  char comments[256];
  sprintf(comments, "%s patched", srcFile);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
