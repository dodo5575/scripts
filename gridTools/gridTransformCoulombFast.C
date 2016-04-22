#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

#include "useful.H"
#include "Grid.H"
#include "BaseGrid.H"
#include "CellDecomposition.H"
#include "Scatter.H"

int countChargeFile(const char* fileName) {
  int nRead;
  int n = 0;
  double q, x, y, z;
  char line[256];

  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("countChargeFile Couldn't open file %s\n.",fileName);
    exit(-1);
  }

  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead = sscanf(line, "%lf %lf %lf %lf", &q, &x, &y, &z);
    if (nRead >= 4) n++;
  }
    
  fclose(inp);
  return n;
}

// Read coordinates into a Vector array.
void readChargeFile(const char* fileName, int num, double* charge, Vector3* pos) {
  int nRead;
  int n = 0;
  double q, x, y, z;
  char line[256];

  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("readChargeFile Couldn't open file %s\n.",fileName);
    exit(-1);
  }

  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead = sscanf(line, "%lf %lf %lf %lf", &q, &x, &y, &z);
    if (nRead >= 4) {
      charge[n] = q;
      pos[n].x = x;
      pos[n].y = y;
      pos[n].z = z;
      n++;
    }
  }
    
  fclose(inp);
}

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
    printf("countListFile Couldn't open file %s\n.",fileName);
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
    printf("readListFile Couldn't open file %s\n.",fileName);
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

// Add Coulomb charges to "dest" using a spatial decomposition.
void addCoulombCoarse(double* chargeQ, Vector3* chargeR, int chargeN, double cutoff, Grid* dest) {
  // Generate the cell decomposition.
  Matrix3 sys = dest->getBox();
  Vector3 origin = dest->getOrigin();
  printf("\nGenerating cell decomposition with cutoff = %.10g...\n", cutoff);
  CellDecomposition cell(sys, origin, cutoff);

  // Add the charges to the cellDecomposition.
  cell.decompose(chargeR, chargeN);
  printf("Cell decomposition of %i cells\n", cell.length());

  // Loop through each cell and get the cell's charge and center of charge.
  // These will serve as the coarse level of electrostatics.
  int cellN = cell.length();
  Vector3* cellR = new Vector3[cellN];
  double* cellQ = new double[cellN];
  for (int c = 0; c < cellN; c++) {
    IndexList local = cell.getCellContents(c);
    
    //Get the center of charge, if it exists.
    double qSum = 0.0;
    Vector3 qrSum = 0.0;
    for (int i = 0; i < local.length(); i++) {
      int ind = local.get(i);
      qSum += chargeQ[ind];
      qrSum += chargeR[ind]*chargeQ[ind];
    }

    // Set the center of charge.
    if (qSum == 0.0) cellR[c] = cell.getPosition(c);
    else cellR[c] = qrSum/qSum;
    // Set the charge of the cell.
    cellQ[c] = qSum;
  }
  printf("Got the centers of charge for %d cells.\n", cellN);

  // Loop through each point of the grid and add Coulomb potentials.
  int n = dest->length();
  int complete = 0;
  for (int i = 0; i < n; i++) {
    // Get the indices (in pos) of possible neighboring points.
    Vector3 r = dest->getPosition(i);
    IndexList neigh = cell.neighborhood(r);
    const int nNeighs = neigh.length();
    
    // Add the local field contributions explicitly.
    for (int j = 0; j < nNeighs; j++) {
      int ind = neigh.get(j);
       
      dest->addCoulomb(chargeQ[ind], chargeR[ind]);
    }

    // Add the distant field contributions with coarse representations.
    for (int c = 0; c < cellN; c++) {
      // Skip neighbors, which we have done explicitly.
      if (neigh.find(c) >= 0) continue;

      dest->addCoulomb(cellQ[c], cellR[c]);
    } 

    // Write the progress.
    int comp = (100*i)/n;
    if (abs(complete - comp) >= 5) {
      printf("%d percent complete\n", comp);
      complete = comp;
    }
  }

  delete[] cellR;
  delete[] cellQ;
}

///////////////////////////////////////////////////////////////////////
// Driver

int mainAdd(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("Add one grid to another applying the transformation given in gridListFile.\n");
    printf("Usage: %s gridListFile coulombConstant cutoff outFile\n", argv[0]);
    printf("gridListFile format: 'gridFile coulombFile exx exy exz eyx eyy eyz ezx ezy ezz ox oy oz'\n");
    printf("coulombFile format: 'q x y z'\n");
    return 0;
  }

  // Count the number of grids.
  int n = countListFile(argv[1]);
  const double coulombConst = strtod(argv[2], NULL);
  double cutoff = strtod(argv[3], NULL);
  printf("Adding %d grids.\n", n);

  // Allocate the memory.
  String* name = new String[n];
  String* coulombFile = new String[n];
  Matrix3* basis = new Matrix3[n];
  Vector3* origin = new Vector3[n];

  // Load the grid names and transformations.
  printf("Reading the list of grid names and transformations from %s.\n", argv[1]);
  readListFile(argv[1], n, name, coulombFile, basis, origin);

  //////// Grids
  // Find the unique grids.
  printf("Finding the unique grids.\n");
  String* uniqueName = new String[n];
  int* gridIndex = new int[n];
  int nGrids = 0;
  for (int i = 0; i < n; i++) {
    // Check if the grid has already been found.
    int alreadyIn = false;
    for (int j = 0; j < nGrids; j++) {
      if (uniqueName[j] == name[i]) {
	alreadyIn = true;
	gridIndex[i] = j;
	break;
      }
    }

    // Add the grid.
    if (!alreadyIn) {
      uniqueName[nGrids] = name[i];
      gridIndex[i] = nGrids;
      nGrids++;
    }
  }

  // Load the unique grids.
  printf("Loading %d unique grids.\n", nGrids);
  Grid** gridList = new Grid*[nGrids];
  for (int i = 0; i < nGrids; i++) {
    printf("Loading %s.\n", uniqueName[i].val());
    gridList[i] = new Grid(uniqueName[i].val());
  }

  // Count the total number of charges.
  int* chargeCount = new int[n];
  int chargeTotal = 0;
  for (int i = 1; i < n; i++) {
    chargeCount[i] = countChargeFile(coulombFile[i]);
    chargeTotal += chargeCount[i];
  }
  printf("Read %d Coulomb charges.\n", chargeTotal);
  double* chargeQ = new double[chargeTotal];
  Vector3* chargeR = new Vector3[chargeTotal];
  int chargeInd = 0;

  // Load the first grid and transform.
  Grid g0(*gridList[gridIndex[0]]);
  Matrix3 b = g0.getBasis();
  Vector3 o = g0.getOrigin();
  g0.setOrigin(basis[0].transform(o) + origin[0]);
  g0.setBasis(basis[0].transform(b));
  printf("Grid %d: %s\n", 0, name[0].val());

  // Interpolate the others onto it.
  for (int i = 1; i < n; i++) {
    printf("Grid %d: %s\n", i, name[i].val());
   
    // Transform the grid.
    Grid* g = gridList[gridIndex[i]];
    Matrix3 b = g->getBasis();
    Vector3 o = g->getOrigin();
    g->setOrigin(basis[i].transform(o) + origin[i]);
    g->setBasis(basis[i].transform(b));

    // Add the grid.
    g0.addInterpolate(*g);

    // Add the coulomb charges to the charge list.
    int points = chargeCount[i];
    double* charge = new double[points];
    Vector3* point = new Vector3[points];
    readChargeFile(coulombFile[i], points, charge, point);
    printf("Registering %d Coulomb charges.\n", points);
    for (int p = 0; p < points; p++) {
      Vector3 pos = basis[i].transform(point[p]) + origin[i];
      
      //g0.addCoulomb(charge[p]*coulombConst, pos);
      chargeQ[chargeInd] = charge[p]*coulombConst;
      chargeR[chargeInd] = pos;
      chargeInd++;
    }
    delete[] charge;
    delete[] point;
        
    // Return the grid to its original state.
    g->setOrigin(o);
    g->setBasis(b);
  }
  
  // Add the Coulomb charges.
  printf("Adding %d Coulomb charges using a spatial decomposition.\n", chargeTotal);
  addCoulombCoarse(chargeQ, chargeR, chargeTotal, cutoff, &g0);

  // Write the result.
  char comments[256];
  snprintf(comments, 256, "%s coulombConst %g cutoff %g", argv[1], coulombConst, cutoff);
  printf("%s\n", comments);
  g0.write(argv[argc-1], comments);

  delete[] name;
  delete[] coulombFile;
  delete[] basis;
  delete[] origin;
  delete[] uniqueName;
  delete[] gridIndex;
  delete[] chargeCount;
  delete[] chargeQ;
  delete[] chargeR;
  for (int i = 0; i < nGrids; i++) delete gridList[i];
  delete[] gridList;
  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
