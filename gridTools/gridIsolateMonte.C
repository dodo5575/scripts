#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
using namespace std;

#include "useful.H"
#include "Grid.H"
#include "RandomGsl.H"

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


void twiddle(Grid* trial, Grid* mask, double jump, Random* randoGen) {
  int n = trial->length();
  
  // Make the Monte Carlo move.
  //for (int i = 0; i < n; i++) {}
  int j = int(floor(n*randoGen->uniform()));
  while (mask->getValue(j) < 0.5) {
    j = int(floor(n*randoGen->uniform()));
  }
      
  double v = trial->getValue(j);
  trial->setValue(j, jump*(randoGen->uniform()-0.5) + v);
}

double error(const Grid* master, const Grid* current, const Grid* mask, const Matrix3* basis, const Matrix3* basisInv, const Vector3* origin, int gridNum) {
  const int n = current->length();

  double error = 0.0;

#pragma omp parallel for reduction(+:error) 
  for (int i = 0; i < n; i++) {
    if (mask->getValue(i) < 0.5) continue;

    double sumVal = current->getValue(i);
    Vector3 r = current->getPosition(i);

    // Transform into the world.
    Vector3 worldR = basis[0].transform(r) + origin[0];

    // Get the master value at this position.
    double masterVal = master->interpolatePotential(worldR);

    // Get the grid indices of the neighboring instances at this position.
    for (int neigh = 1; neigh < gridNum; neigh++) {
      // Convert the position into neighbors frame.
      Vector3 neighR = basisInv[neigh].transform(worldR - origin[neigh]);
      
      sumVal += current->interpolatePotential(neighR);
    }

    // Add the error.
    error += (masterVal-sumVal)*(masterVal-sumVal);
  }

  return error;
}

///////////////////////////////////////////////////////////////////////
// Driver

int mainAdd(int argc, char* argv[]) {
  if ( argc != 10 ) {
    printf("Add one grid to another applying the transformation given in gridListFile.\n");
    printf("Usage: %s masterGrid gridListFile initGrid maskGrid monteSteps temperature jump outFile\n", argv[0]);
    printf("gridListFile format: 'gridFile coulombFile exx exy exz eyx eyy eyz ezx ezy ezz ox oy oz'\n");
    return 0;
  }

  // Load the master grid.
  Grid master(argv[1]);
  const char* listFile = argv[2];
  const char* initGrid = argv[3];
  const char* maskGrid = argv[4];
  const int monteSteps = atoi(argv[5]);
  const int checkSteps = atoi(argv[6]);
  const double temper0 = strtod(argv[7], NULL);
  const double jump0 = strtod(argv[8], NULL);
  const char* outFile = argv[argc-1];

  char previewFile[256];
  snprintf(previewFile, 256, "preview_%s", outFile);

  // Make the random number generator.
  long unsigned int seed = (unsigned int)time((time_t *)NULL);
  Random randoGen(seed);

  // Count the number of grids.
  int n = countListFile(listFile);
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
  printf("Reading the list of grid names and transformations from %s.\n", listFile);
  readListFile(listFile, n, name, coulombFile, basis, origin);
  for (int i = 0; i < n; i++) basisInv[i] = basis[i].inverse();

  // Open the first grid.
  Grid init(initGrid);
  Grid mask(maskGrid);
  //const int initSize = init.length();

  // Perform Monte Carlo steps.
  int accept = 0;
  double jump = jump0;
  double e0 = error(&master, &init, &mask, basis, basisInv, origin, n);
  double bestError = e0;
  Grid best(init);
  bool bestUpdate = false;
  for (int m = 1; m <= monteSteps; m++) {
    // Simulated annealing.
    double temper = temper0*(1.0-double(m-1)/monteSteps);
    
    // Make the trial move.
    Grid trial(init);
    twiddle(&trial, &mask, jump, &randoGen);

    // Check the energy.
    double e1 = error(&master, &trial, &mask, basis, basisInv, origin, n);
    
    // Metropolis
    if (e1 < e0 || exp((e0-e1)/temper) > randoGen.uniform()) {
      // Accept the move.
      init = trial;
      accept++;
      e0 = e1;

      // Keep track of the best configuration.
      if (e0 < bestError) {
	bestError = e0;
	best = init;
	//printf("BEST: step %d error %g\n", m, bestError);
	bestUpdate = true;
      }
    }
      
    // Check.
    if (m % checkSteps == 0) {
      double ratio = double(accept)/checkSteps;
      accept = 0;
      
      printf("step %d error %g jump %g temper %g acceptance %g\n", m, e0, jump, temper, ratio);

      if (ratio < 0.5) jump /= 1.0+randoGen.uniform();
      else jump *= 1.0+randoGen.uniform();

      // Write the best configuration.
      if (bestUpdate) {
	bestUpdate = false;
	best.write(previewFile);
      }
    }
  }

  // Done with the Monte Carlo steps.
  printf("Done.\n");
  best.write(outFile);
    
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
