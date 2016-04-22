#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

#include "useful.H"
#include "Grid.H"

int countListFile(const char* fileName) {
  int nRead;
  int n = 0;
  double exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz, ox, oy, oz;
  char line[1024];
  char word0[1024];
  char word1[1024];
  //char word2[256];

  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
    exit(-1);
  }

  while (fgets(line, 1024, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead = sscanf(line, "%s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		   word0, word1, &exx, &exy, &exz, &eyx, &eyy, &eyz, &ezx, &ezy, &ezz, &ox, &oy, &oz);
    if (nRead >= 14) n++;
  }
    
  fclose(inp);
  return n;
}


// Read coordinates into a Vector array.
void readListFile(const char* fileName, int num, String* name0, String* name1, Matrix3* basis, Vector3* origin) {
 int nRead;
  int n = 0;
  double exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz, ox, oy, oz;
  char line[1024];
  char word0[1024];
  char word1[1024];
  //char word2[256];

  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
    exit(-1);
  }

  while (fgets(line, 1024, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead = sscanf(line, "%s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		   word0, word1, &exx, &exy, &exz, &eyx, &eyy, &eyz, &ezx, &ezy, &ezz, &ox, &oy, &oz);
    if (nRead >= 13) {
      name0[n] = word0;
      name1[n] = word1;
      //name2[n] = word2;
      basis[n].exx = exx;
      basis[n].exy = exy;
      basis[n].exz = exz;
      basis[n].eyx = eyx;
      basis[n].eyy = eyy;
      basis[n].eyz = eyz;
      basis[n].ezx = ezx;
      basis[n].ezy = ezy;
      basis[n].ezz = ezz;
      origin[n].x = ox;
      origin[n].y = oy;
      origin[n].z = oz;
      n++;
    }
  }
    
  fclose(inp);
}


///////////////////////////////////////////////////////////////////////
// Driver

int mainAdd(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("Add one grid to another applying the transformation given in gridListFile.\n");
    printf("Usage: %s masterGrid gridListFile iterations outFile\n", argv[0]);
    printf("gridListFile format: 'gridFile coulombFile exx exy exz eyx eyy eyz ezx ezy ezz ox oy oz'\n");
    return 0;
  }

  // Load the master grid.
  Grid master(argv[1]);
  int iterNum = atoi(argv[3]);

  // Count the number of grids.
  int n = countListFile(argv[2]);
  printf("Found %d grids.\n", n);
  int neighborNum = n-1;
  printf("Found %d neighbors.\n", neighborNum);
  
  // Allocate the memory.
  String* name = new String[n];
  String* coulombFile = new String[n];
  Matrix3* basis = new Matrix3[n];
  Matrix3* basisInv = new Matrix3[n];
  Vector3* origin = new Vector3[n];

  // Load the grid names and transformations.
  printf("Reading the list of grid names and transformations from %s.\n", argv[1]);
  readListFile(argv[2], n, name, coulombFile, basis, origin);
  for (int i = 0; i < n; i++) basisInv[i] = basis[i].inverse();

  // Initial conditions.
  printf("Initial conditions from `%s'\n", name[0].val());
  Grid init(name[0].val());
  const int initSize = init.length();
  Grid iter(init);

  Grid* last = &init;
  Grid* next = &iter;

  Grid attempt(master);
  attempt.zero();

  for (int i = 0; i < attempt.length(); i++) {
    Vector3 worldR = attempt.getPosition(i);

    for (int neigh = 0; neigh < n; neigh++) {
      // Transform into local space.
      Vector3 localR = basisInv[neigh].transform(worldR - origin[neigh]);
      if (init.inGrid(localR)) {
	double v = init.interpolatePotential(localR);
	attempt.setValue(i,v + attempt.getValue(i));
      }
    }
  }
  char outName[256];
  snprintf(outName, 256, "attempt_%s", argv[argc-1]);
  attempt.write(outName);

  printf("Running %d iterations.\n", iterNum);
  for (int i = 0; i < iterNum; i++) {
    printf("iteration %d\n", i);
    if (i % 1 == 0) {
      char outName[256];
      snprintf(outName, 256, "iter%d_%s", i, argv[argc-1]);
      last->write(outName);
    }

    // Loop through the component grid.
    for (int j = 0; j < initSize; j++) {
      Vector3 r = init.getPosition(j);
      // Transform into the world.
      Vector3 worldR = basis[0].transform(r) + origin[0];

      // Get the master value at this position.
      double val = master.interpolatePotential(worldR);

      // Now subtract the neighbor values.
      for (int neigh = 1; neigh < n; neigh++) {
	// Convert the position into neighbors frame.
	Vector3 neighR = basisInv[neigh].transform(worldR - origin[neigh]);
	// Subtract the value.
	if (init.inGrid(neighR)) val -= last->interpolatePotential(neighR);
      }

      // Make the new grid.
      next->setValue(j, val);
    }

    // Swap the pointers.
    Grid* tmp = last;
    last = next;
    next = tmp;
  }

  // Write the result.
  char comments[256];
  snprintf(comments, 256, "%s sum", argv[2]);
  printf("%s\n", comments);
  last->write(argv[argc-1], comments);

  delete[] name;
  delete[] coulombFile;
  delete[] basis;
  delete[] origin;
  delete[] basisInv;
  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
