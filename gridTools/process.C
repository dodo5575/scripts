#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../useful.H"
#include "./Grid.H"

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

  if ( argc != 3 ) {
    printf("Usage: %s inFile outFile\n", argv[0]);
    return 0;
  }

  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  // Write the profile.
  sprintf(s, "%s", argv[2]);
  printf("Writing the profile to %s...\n", s);
  g.averageProfileZBoltzmann(s);

  return 0;
}

//Rogan
int mainProfileAverage(int argc, char* argv[])
{
  char s[128];

  if ( argc != 3 ) {
    printf("Usage: %s inFile outFile\n", argv[0]);
    return 0;
  }

  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  // Write the profile.
  sprintf(s, "%s", argv[2]);
  printf("Writing the profile to %s...\n", s);
  g.averageProfileZ(s);

  return 0;
}

int mainProfileAverageDeviation(int argc, char* argv[])
{
  char s[128];

  if ( argc != 3 ) {
    printf("Usage: %s inFile outFile\n", argv[0]);
    return 0;
  }

  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  // Write the profile.
  sprintf(s, "%s", argv[2]);
  printf("Writing the profile to %s...\n", s);
  g.averagedeviationProfileZ(s);

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

int mainAverageRegion(int argc, char* argv[])
{
  if ( argc != 5 ) {
    printf("Usage: %s inFile outFile z0 z1\n", argv[0]);
    return 0;
  }

  double z0 = strtod(argv[3], NULL);
  double z1 = strtod(argv[4], NULL);

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

  char comments[256];
  sprintf(comments, "%s shifted by %.10g", argv[1], -v);
  printf("%s\n", comments);
  g.write(argv[2], comments);

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

  return 1;
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

  return 1;
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
    //double fz = -1;

    fprintf(out, "%.10g %.10g\n", z, f.z);
  }
  fclose(out);

  return 1;
}

int mainShiftToAverage(int argc, char* argv[])
{
  char s[128];

  if ( argc != 7 ) {
    printf("Usage: %s inFile outFile x y zmin zmax\n", argv[0]);
    return 0;
  }

  double x = strtod(argv[3], NULL);
  double y = strtod(argv[4], NULL);
  double zA = strtod(argv[5], NULL);
  //double zB = strtod(argv[6], NULL);

  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  Vector3 r0 = g.getOrigin();
  Matrix3 b = g.getBasis();

  int ix = int ((x - r0.x) / b.exx);
  int iy = int ((y - r0.y) / b.eyy);
  int z0 = int ((zA - r0.z) / b.ezz);
  int z1 = int ((zA - r0.z) / b.ezz);

  printf("Averaging along\n\tx = %.2f\ty = %.2f\nfrom\n\tz0 = %.2f\tto\nz0 = %.2f\n", ix*b.exx+r0.x, iy*b.eyy+r0.y, z0*b.ezz+r0.z, z1*b.ezz+r0.z);

  // Shift the average over the top of the grid to zero.
  double v0 = g.averageZLine(ix, iy, z0, z1);
  g.shift(-v0);
  printf("Shifting by %.10g\n", -v0);

  // Write the profile down the center and the boltzmann average profile.
  sprintf(s, "%s.z.dat", argv[2]);
  g.profileZ(0.5, 0.5, s);
  sprintf(s, "%s.boltzz.dat", argv[2]);
  g.averageProfileZBoltzmann(s);

  // Write the result.
  g.write(argv[2], "");
  printf("\nWrote %s.\n", argv[2]);
  return 0;
}


int mainShiftToAverageOfBox(int argc, char* argv[])
{
  char s[128];

  if ( argc != 7 ) {
    printf("Usage: %s inFile outFile x y zmin zmax\n", argv[0]);
    return 0;
  }

  double x = strtod(argv[3], NULL);
  double y = strtod(argv[4], NULL);
  double zA = strtod(argv[5], NULL);
  double zB = strtod(argv[6], NULL);

  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  Vector3 r0 = g.getOrigin();
  Matrix3 b = g.getBasis();

  int ix = int ((x - r0.x) / b.exx);
  int iy = int ((y - r0.y) / b.eyy);
  int z0 = int ((zA - r0.z) / b.ezz);
  int z1 = int ((zB - r0.z) / b.ezz);

  printf("Averaging \n\tx = %.2f, %.2f\ty = %.2f, %.2f\n\tz0 = %.2f, %.2f\n", (ix-1)*b.exx+r0.x, (ix+1)*b.exx+r0.x, (iy-1)*b.eyy+r0.y, (iy+1)*b.eyy+r0.y, z0*b.ezz+r0.z, z1*b.ezz+r0.z);

  // Shift the average over the top of the grid to zero.
  double v0 = g.averageZBox(ix-1, ix+1, iy-1, iy+1, z0, z1);
  g.shift(-v0);
  printf("Shifting by %.10g\n", -v0);

  // Write the profile down the center and the boltzmann average profile.
  sprintf(s, "%s.z.dat", argv[2]);
  g.profileZ(0.5, 0.5, s);
  sprintf(s, "%s.boltzz.dat", argv[2]);
  g.averageProfileZBoltzmann(s);

  // Write the result.
  g.write(argv[2], "");
  printf("\nWrote %s.\n", argv[2]);
  return 0;
}

int mainShiftToAverageExact(int argc, char* argv[])
{
  char s[128];

  if ( argc != 7 ) {
    printf("Usage: %s inFile outFile x y zmin zmax\n", argv[0]);
    return 0;
  }

  double x = strtod(argv[3], NULL);
  double y = strtod(argv[4], NULL);
  double zA = strtod(argv[5], NULL);
  double zB = strtod(argv[6], NULL);

  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  Vector3 r0 = g.getOrigin();
  Matrix3 b = g.getBasis();

  int ix = int ((x - r0.x) / b.exx);
  double dx1 = (x-(ix*b.exx+r0.x))/b.exx;
  double dx2 = 1-dx1;

  int iy = int ((y - r0.y) / b.eyy);
  double dy1 = (y-(iy*b.eyy+r0.y))/b.eyy;
  double dy2 = 1-dy1;

  //  printf("X: x %.2f, x0  %.2f, x1 %.2f, dx0 %.2f, dx1, %.2f\n", x, ix*b.exx+r0.x, (ix+1)*b.exx+r0.x, xD0, xD1);
  //  printf("Y: y %.2f, y0  %.2f, y1 %.2f, dy0 %.2f, dy1, %.2f\n", y, iy*b.eyy+r0.y, (iy+1)*b.eyy+r0.y, yD0, yD1);

  int z0 = int ((zA - r0.z) / b.ezz);
  int z1 = int ((zB - r0.z) / b.ezz);

  printf("Averaging along\tx = %.2f\ty = %.2f\nfrom\t\tz0 = %.2f\tto\tz0 = %.2f\n", ix*b.exx+r0.x, iy*b.eyy+r0.y, z0*b.ezz+r0.z, z1*b.ezz+r0.z);


  double errorx1, errorx2, errory1, errory2;
  double stddevx1, stddevx2, stddevy1, stddevy2;
  /* Bilinear mix */
  // Mix X
  double vya = dx2*g.averageZLine(ix, iy, z0, z1, errorx1, stddevx1) + dx1*g.averageZLine(ix+1, iy, z0, z1, errorx2, stddevx2);
  double vyb = dx2*g.averageZLine(ix, iy+1, z0, z1, errory1, stddevy1) + dx1*g.averageZLine(ix+1, iy+1, z0, z1, errory2, stddevy2);
  // Mix Y
  double v0 = dy2*vya + dy1*vyb;

  double error = dy2*(dx2*errorx1+dx1*errorx2)+dy1*(dx2*errory1+dx1*errory2);
  double stddev = dy2*(dx2*stddevx1+dx1*stddevx2)+dy1*(dx2*stddevy1+dx1*stddevy2);

  g.shift(-v0);
  printf("Shifting by %.10g.  Error in graph is %.10g / %.10g.\n", -v0, error, stddev);

  // Write the profile down the center and the boltzmann average profile.
  sprintf(s, "%s.z.dat", argv[2]);
  g.profileZ(0.5, 0.5, s);
  sprintf(s, "%s.boltzz.dat", argv[2]);
  g.averageProfileZBoltzmann(s);

  // Write the result.
  g.write(argv[2], "");
  printf("\nWrote %s.\n", argv[2]);
  return 0;
}

int mainPlotZ(int argc, char* argv[])
{
  //char s[128];

  if ( argc != 7 ) {
    printf("Usage: %s inFile outFile x y zmin zmax\n", argv[0]);
    return 0;
  }

  double x = strtod(argv[3], NULL);
  double y = strtod(argv[4], NULL);
  double z0 = strtod(argv[5], NULL);
  double z1 = strtod(argv[6], NULL);

  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  g.profileZ(x,y,z0,z1,argv[2]);

  printf("\nWrote %s.\n", argv[2]);
  return 0;
}

// Just crop along the Z
int mainCropZ(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("Usage: %s inFile outFile z0 z1\n", argv[0]);
    return 0;
  }

  // Load the first grid.
  Grid g(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  int z0 = g.nearestZIndexUp(strtod(argv[3], NULL));
  int z1 = g.nearestZIndexDown(strtod(argv[4], NULL));
  g.crop(z0,z1);

  // Write the result.
  printf("Writing the cropped grid to %s.\n", argv[2]);
  g.write(argv[2]);

  return 1;
}

int main(int argc, char* argv[]) {
  //return mainProfileBoltzmann(argc, argv);
  //return mainProfileAverage(argc,argv); //Rogan
  //return mainRoganNegativeAverage(argc,argv); //Rogan
  //return mainProfileX(argc, argv);
  //return mainShiftToAverage(argc, argv);
  //return mainShiftToAverageOfBox(argc, argv);
  //return mainShiftToAverageExact(argc, argv);
  //mainCrop(argc, argv);
  //return mainPlotZ(argc,argv);
  //return mainAverageRegion(argc,argv);
  //return mainProfileAverageDeviation(argc,argv);
  return mainCropZ(argc, argv);
}
