#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

#include "useful.H"
#include "Grid.H"

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
  if ( argc != 4 ) {
    printf("Extract the data from the single column data file and insert it in a grid with dimensions of inGrid.\n");
    printf("Usage: %s inGrid dataFile outGrid\n", argv[0]);
    return 0;
  }

  const char* inGrid = argv[1];
  const char* dataFile = argv[2];
  const char* outGrid = argv[argc-1];

  // Load the grid.
  Grid src(inGrid);
  const int size = src.length();
  printf("Loaded `%s' with %d nodes.\n", inGrid, size);


  // Open the file.
  FILE* inp = fopen(dataFile,"r");
  if (inp == NULL) {
    printf("ERROR: Couldn't open file `%s'.\n", dataFile);
    exit(-1);
  }

  // Insert its data.
  char line[256];
  double v;
  int node = 0;
  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 1) continue;

    // Put this value in the grid.
    int nRead = sscanf(line, "%lf", &v);
    if (nRead == 1) {
      if (node < size) src.setValue(node, v);
      node++;
    }
  }

  // Check number read.
  if (node > size) 
    printf("Warning: Too many data values! Read %d for %d grid points.\n", node, size);
  if (node < size)
    printf("Warning: Too few data values! Read %d for %d grid points.\n", node, size);

  // Write the result.
  char comments[256];
  snprintf(comments, 256, "%s with data from %s", inGrid, dataFile);
  printf("%s\n", comments);
  src.write(outGrid, comments);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
