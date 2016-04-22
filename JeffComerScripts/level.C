#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver
void add(Grid& g, const char* fileName) {
  Grid g1(fileName);
  
  g.addGrid(g1);
}

int mainAverage(int argc, char* argv[])
{
  char s[128];
  // Factor to obtain volts from e/kT
  double factor = 0.0258520296;
  //double factor = 1.0;

  if ( argc < 4 ) {
    printf("Usage: %s outFile voltage dxFile0 dxFile1 ... \n", argv[0]);
    return 0;
  }
  double voltage = strtod(argv[2],NULL);
  
  int n = argc-3;
  printf("\nPROCESSING GRIDS\n");
  printf("%d grids will be averaged\n", n);
  printf("Scale factor: %.10g\n", factor);
  printf("Voltage drop: %.10g\n", voltage);
  printf("Output file: %s\n", argv[1]);
  

  // Load the first grid.
  Grid g(argv[3]);
  printf("\nLoaded %s.\n", argv[3]);

  printf("Adding grids...\n");
  for (int i = 4; i < argc; i++) {
    add(g, argv[i]);
    printf("Loaded %s.\n", argv[i]);
  }

  // Scale to make an average and convert the units.
  //printf("\nConverting from kT/e to volts...\n");
  g.scale(factor/(double)n);

  // Compute the external field.
  printf("\nApplying the external field...\n");
  Vector3 extent = g.getExtent();  
  Vector3 eField(0.0, 0.0, -voltage/extent.z);
  printf("System length along z:  %.10g\n", extent.z);
  printf("voltage = %.10g => eField = %s\n", voltage, (const char*)(eField.toString()));
  g.addGradient(eField);

  // Shift the average over the top of the grid to zero.
  double v0 = g.averageSection(0);
  g.shift(-v0);
  printf("Shifting by %.10g\n", -v0);

  // Write the profile down the center and the average profile.
  sprintf(s, "%s.z.dat", argv[1]);
  g.profileZ(0.5, 0.5, s);
  sprintf(s, "%s.avgz.dat", argv[1]);
  g.averageProfileZ(s);

  // Write the result.
  g.write(argv[1], "");
  printf("\nWrote %s.\n", argv[1]);
  return 0;
}

int mainProfile(int argc, char* argv[])
{
  char s[128];

  if ( argc != 3 ) {
    printf("Usage: %s inFile factor\n", argv[0]);
    return 0;
  }
  
  double factor = strtod(argv[2], NULL);

  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  // Write the profile.
  sprintf(s, "%s.dat", argv[1]);
  printf("Writing the profile to %s...\n", s);
  g.profileZ(factor, factor, s);
  
  return 0;
}

int mainProfileBoltzmann(int argc, char* argv[])
{
  char s[128];

  if ( argc != 2 ) {
    printf("Usage: %s inFile\n", argv[0]);
    return 0;
  }
  
  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  // Write the profile.
  sprintf(s, "%s.dat", argv[1]);
  printf("Writing the profile to %s...\n", s);
  g.averageProfileZBoltzmann(s);
  
  return 0;
}

int mainProfileCyl(int argc, char* argv[])
{
  //const double factor = 1.660538782/15.9994*18.0154; // (g/cm^3)/(Da/A^3)
  const double factor = 1.0;
  char s[128];

  if ( argc != 3 ) {
    printf("Usage: %s inFile radius\n", argv[0]);
    return 0;
  }
  
  double radius = strtod(argv[2], NULL);

  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  g.scale(factor);
  printf("\nScaled by %.10g.\n", factor);

  // Write the profile.
  sprintf(s, "%s.avg.dat", argv[1]);
  printf("Writing the profile to %s...\n", s);
  g.cylinderProfileZ(s, radius);
  
  return 0;
}

int mainAverageRegion1(int argc, char* argv[]) {
  if ( argc != 8 ) {
    printf("Usage: %s inFile x0 y0 z0 x1 y1 z1\n", argv[0]);
    return 0;
  }
  
  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);
  
  Vector3 r0;
  Vector3 r1;
  
  r0.x = strtod(argv[2],NULL);
  r0.y = strtod(argv[3],NULL);
  r0.z = strtod(argv[4],NULL);
  r1.x = strtod(argv[5],NULL);
  r1.y = strtod(argv[6],NULL);
  r1.z = strtod(argv[7],NULL);

  double v = g.averageRegion(r0,r1);

  g.shift(-v);
  printf("Shifted by %.10g.\n", -v);
  
  char s[256];
  sprintf(s, "shift_%s", argv[1]);
  printf("%s\n", s);
  char comments[256];
  sprintf(comments, "%s shifted by %.10g", argv[1], -v);
  printf("%s\n", comments);
  g.write(s, comments);
  
  return 0;
}

int mainAverageRegion(int argc, char* argv[])
{
  if ( argc != 4 ) {
    printf("Usage: %s inFile z0 z1\n", argv[0]);
    return 0;
  }

  double z0 = strtod(argv[2], NULL);
  double z1 = strtod(argv[3], NULL);

  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  Vector3 r0 = 2.0*g.getOrigin();
  Vector3 r1 = 2.0*(g.getOrigin()+g.getExtent());
  r0.z = z0;
  r1.z = z1;
  double v = g.averageRegion(r0,r1);

  g.shift(-v);
  printf("Shifted by %.10g.\n", -v);
  
  char s[256];
  sprintf(s, "shift_%s", argv[1]);
  printf("%s\n", s);
  char comments[256];
  sprintf(comments, "%s shifted by %.10g", argv[1], -v);
  printf("%s\n", comments);
  g.write(s, comments);
  
  return 0;
}

 int mainDepth(int argc, char* argv[])
{
  if ( argc != 2 ) {
    printf("Usage: %s inFile\n", argv[0]);
    return 0;
  }
 
  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  // Write the depth map.
  char s[256];
  sprintf(s, "%s.depth", argv[1]);
  printf("Writing the depth map to %s.\n", s);
  g.depthMapZ(s);

  return 0;
}

int mainCrop(int argc, char* argv[]) {
  if ( argc != 8 ) {
    printf("Usage: %s inFile x0 y0 z0 x1 y1 z1\n", argv[0]);
    return 0;
  }
  
  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);
  
  Vector3 r0;
  Vector3 r1;
  
  r0.x = strtod(argv[2],NULL);
  r0.y = strtod(argv[3],NULL);
  r0.z = strtod(argv[4],NULL);
  r1.x = strtod(argv[5],NULL);
  r1.y = strtod(argv[6],NULL);
  r1.z = strtod(argv[7],NULL);
  
  printf("Keeping from (%s) to (%s).\n", r0.toString().val(), r1.toString().val());
  g.crop(r0,r1);
  
  // Write the result.
  char comments[256];
  sprintf(comments, "%s cropped to (%s) (%s)", argv[1], r0.toString().val(), r1.toString().val());
  char s[256];
  sprintf(s, "crop_%s", argv[1]);
  printf("Writing the cropped grid to %s.\n", s);
  g.write(s, comments);
}

int mainList(int argc, char* argv[]) {
  if ( argc != 2 ) {
    printf("Usage: %s inFile\n", argv[0]);
    return 0;
  }
  
  // Load the grid
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);
  
  char outName[256];
  sprintf(outName, "%s.dat", argv[1]);
  printf("Writing grid to %s.\n", outName);
  g.writeData(outName);
}

int mainForce(int argc, char* argv[]) {
  if ( argc != 4 ) {
    printf("Usage: %s gridFile z0 z1\n", argv[0]);
    return 0;
  }

  double z0 = strtod(argv[2],NULL);
  double z1 = strtod(argv[3],NULL);
  char outName[256];
  sprintf(outName, "%s.dat", argv[1]); 

  const int n = 2000;
  double dz = (z1 - z0)/n;
  Grid g(argv[1]);
  
  FILE* out = fopen(outName, "w");
  for (int i = 0; i < n; i++) {
    double z = z0 + i*dz;
    //double v = g.interpolatePotential(Vector3(0,0,z));
    Vector3 f = g.interpolateForce(Vector3(0,0,z));
    //double fz = g.interpolateGrad(Vector3(z,0,0), 0);
    double fz = -1;

    fprintf(out, "%.10g %.10g\n", z, f.z);
  }
  fclose(out);
}

int mainCropValid(int argc, char* argv[]) {
  if ( argc != 4 ) {
    printf("Usage: %s cropGrid validGrid outgrid\n", argv[0]);
    return 0;
  }
  
  Grid crop(argv[1]);
  Grid valid(argv[2]);

  //crop.crop(valid);

  char comments[256];
  sprintf(comments, "Cropped %s using %s.", argv[1], argv[2]);
  crop.write(argv[3], comments);
}

int main(int argc, char* argv[]) {
  return mainAverageRegion1(argc, argv);
}
