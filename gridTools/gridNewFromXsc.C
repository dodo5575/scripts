#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"
#include "Scatter.H"

///////////////////////////////////////////////////////////////////////
// Driver

Matrix3 readXsc(const char* xscFile) {
  // Open the file.
  FILE* inp = fopen(xscFile, "r");
  if (inp == NULL) {
    printf("readXsc Couldn't open file %s\n.", xscFile);
    exit(-1);
  }
  
  char line[256];
  int steps;
  double exx, eyx, ezx, exy, eyy, ezy, exz, eyz, ezz;
  while (fgets(line, 256, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
    
      // Read values.
      int nRead = sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &steps, &exx, &eyx, &ezx, &exy, &eyy, &ezy, &exz, &eyz, &ezz);
      printf("%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", steps, exx, eyx, ezx, exy, eyy, ezy, exz, eyz, ezz );
      if (nRead == 10) {
	break;
      }
    }
    fclose(inp);

    return Matrix3(exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz);
}

int main(int argc, char* argv[]) {
  if ( argc != 4 ) {
    printf("Make a new grid (all zeros) from an xsc file and the desired grid spacing.\n");
    printf("\nUsage: %s xscFile gridSpacing outGrid\n", argv[0]);
    return 0;
  }

  const char* xscFile = argv[1];
  double gridSpacing = strtod(argv[2], NULL);
  const char* outGrid = argv[argc-1];
  

  // Get system vectors.
  Matrix3 box = readXsc(xscFile);
  Grid orig(box, gridSpacing);

  char comments[256];
  snprintf(comments, 256, "grid with basis %s ", argv[1]);
  orig.write(outGrid, comments);
  printf("Wrote `%s'.\n", outGrid);

  return 0;
}
