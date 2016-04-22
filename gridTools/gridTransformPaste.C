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
  char word[1024];

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
    nRead = sscanf(line, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		   word, &exx, &exy, &exz, &eyx, &eyy, &eyz, &ezx, &ezy, &ezz, &ox, &oy, &oz);
    if (nRead == 13) n++;
    else printf("countListFile Improperly formatted list file.\n");
  }
    
  fclose(inp);
  return n;
}

// Read coordinates into a Vector array.
void readListFile(const char* fileName, int num, String* name, Matrix3* basis, Vector3* origin) {
 int nRead;
  int n = 0;
  double exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz, ox, oy, oz;
  char line[1024];
  char word[1024];

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
    nRead = sscanf(line, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		   word, &exx, &exy, &exz, &eyx, &eyy, &eyz, &ezx, &ezy, &ezz, &ox, &oy, &oz);
    if (nRead == 13) {
      name[n] = word;
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
  if ( argc != 3 ) {
    printf("Add one grid to another applying the transformation given in gridListFile.\n");
    printf("Usage: %s gridListFile outFile\n", argv[0]);
    printf("gridListFile format: 'gridFileName exx exy exz eyx eyy eyz ezx ezy ezz ox oy oz'\n");
    return 0;
  }

  // Count the number of grids.
  int n = countListFile(argv[1]);
  printf("Adding %d grids.\n", n);

  // Allocate the memory.
  String* name = new String[n];
  Matrix3* basis = new Matrix3[n];
  Vector3* origin = new Vector3[n];

  // Load the grid names and transformations.
  printf("Reading the list of grid names and transformations from %s.\n", argv[1]);
  readListFile(argv[1], n, name, basis, origin);

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
    //g0.addInterpolate(*g);
    g0.pasteInterpolate(*g);

    // Return the grid to its original state.
    g->setOrigin(o);
    g->setBasis(b);
  }

  // Write the result.
  char comments[1024];
  snprintf(comments, 1024, "%s sum", argv[1]);
  printf("%s\n", comments);
  g0.write(argv[argc-1], comments);

  delete[] name;
  delete[] basis;
  delete[] origin;
  delete[] uniqueName;
  delete[] gridIndex;
  for (int i = 0; i < nGrids; i++) delete gridList[i];
  delete[] gridList;
  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
