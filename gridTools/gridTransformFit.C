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
  char word0[256];
  char word1[256];

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
  char word0[256];
  char word1[256];

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
    if (nRead >= 14) {
      name0[n] = word0;
      name1[n] = word1;
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
  if ( argc != 4 ) {
    printf("Paste one grid into another applying the transformation given in gridListFile.\n");
    printf("Usage: %s gridListFile fitThreshold outFile\n", argv[0]);
    printf("gridListFile format: 'gridFileName maskFileName exx exy exz eyx eyy eyz ezx ezy ezz ox oy oz'\n");
    return 0;
  }

  const int blurLevel = 7;
  double fitThreshold = strtod(argv[2], NULL);

  // Count the number of grids.
  int n = countListFile(argv[1]);
  printf("Adding %d grids.\n", n);

  // Allocate the memory.
  String* name = new String[n];
  String* maskName = new String[n];
  Matrix3* basis = new Matrix3[n];
  Vector3* origin = new Vector3[n];

  // Load the grid names and transformations.
  printf("Reading the list of grid names and transformations from %s.\n", argv[1]);
  readListFile(argv[1], n, name, maskName, basis, origin);

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
  
//////// Masks
// Find the unique masks.
  printf("Finding the unique masks.\n");
  String* maskUniqueName = new String[n];
  int* maskIndex = new int[n];
  int nMasks = 0;
  for (int i = 0; i < n; i++) {
    // Check if the grid has already been found.
    int alreadyIn = false;
    for (int j = 0; j < nMasks; j++) {
      if (maskUniqueName[j] == maskName[i]) {
	alreadyIn = true;
	maskIndex[i] = j;
	break;
      }
    }

    // Add the mask.
    if (!alreadyIn) {
      maskUniqueName[nMasks] = maskName[i];
      maskIndex[i] = nMasks;
      nMasks++;
    }
  }

  // Load the unique masks.
  printf("Loading %d unique masks.\n", nMasks);
  Grid** maskList = new Grid*[nMasks];
  for (int i = 0; i < nMasks; i++) {
    printf("Loading %s.\n", maskUniqueName[i].val());
    maskList[i] = new Grid(maskUniqueName[i].val());
  }

  // Load the first grid and transform.
  Grid g0 (*gridList[gridIndex[0]]);
  //Grid m0(*maskList[gridIndex[0]]); Doesn't do anything.
  Matrix3 b = g0.getBasis();
  Vector3 o = g0.getOrigin();
  g0.setOrigin(basis[0].transform(o) + origin[0]);
  g0.setBasis(basis[0].transform(b));
  printf("Grid %d: %s\n", 0, name[0].val());
  Grid result(g0);

  // Interpolate the others onto it.
  for (int i = 1; i < n; i++) {
    printf("Grid %d: %s\n", i, name[i].val());
   
    // Transform the grid.
    Grid g(*gridList[gridIndex[i]]);
    Matrix3 gb = g.getBasis();
    Vector3 go = g.getOrigin();
    g.setOrigin(basis[i].transform(go) + origin[i]);
    g.setBasis(basis[i].transform(gb));
    
    // Transform the mask.
    Grid m(*maskList[maskIndex[i]]);
    Matrix3 mb = m.getBasis();
    Vector3 mo = m.getOrigin();
    m.setOrigin(basis[i].transform(mo) + origin[i]);
    m.setBasis(basis[i].transform(mb));
    
    // Make a boundary mask.
    Grid bord(m);
    bord.alphaThreshold(0.6);
    for (int i = 0; i < blurLevel; i++) bord.blur();
    Grid inv(bord);
    inv.alphaInvert();
    bord.multiply(inv);
    bord.alphaThreshold(0.0001);
    bord.multiply(m);
    bord.alphaThreshold(0.0001);

    // Remove points from the boundary mask that have higher values than the fitThreshold in either grid.
    Grid fitGrid(g);
    fitGrid.alphaThreshold(fitThreshold);
    fitGrid.alphaInvert();
    //fitGrid.write("fitGrid.dx");
    Grid fitGrid0(g0);
    fitGrid0.alphaThreshold(fitThreshold);
    fitGrid0.alphaInvert();
    //fitGrid0.write("fitGrid0.dx");

    bord.multiplyInterpNoWrap(fitGrid);
    bord.multiplyInterpNoWrap(fitGrid0);
    bord.write("boundary.dx");

    // Match the grid on the boundary between the two grids.
    double v0 = g0.averageRegion(bord);
    double v = g.averageRegion(bord);
    g.shift(v0-v);
    
    // Add the grid.
    result.paste(g, m);
  }

  // Write the result.
  char comments[256];
  sprintf(comments, "%s sum", argv[1]);
  printf("%s\n", comments);
  result.write(argv[argc-1], comments);

  delete[] name;
  delete[] maskName;
  delete[] basis;
  delete[] origin;
  delete[] uniqueName;
  delete[] maskUniqueName;
  delete[] gridIndex;
  delete[] maskIndex;
  for (int i = 0; i < nGrids; i++) delete gridList[i];
  delete[] gridList;
  for (int i = 0; i < nMasks; i++) delete maskList[i];
  delete[] maskList;
  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
