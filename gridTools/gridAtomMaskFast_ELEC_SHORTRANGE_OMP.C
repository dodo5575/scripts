// Author: Jeff Comer <jcomer2@illinois.edu>
// Modified for LJ potential: David Wells <dbwells2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cctype>

#include "useful.H"
#include "BaseGrid.H"
#include "Grid.H"
#include "CellDecomposition.H"
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

bool within(const Scatter& pos, Vector3 r, double rMax, const Grid& g) {
  const int n = pos.length();
  double r2 = rMax*rMax;
  for (int p = 0; p < n; p++) {
    Vector3 d = g.wrapDiff(pos.get(p) - r);
    if (d.length2() < r2) return true;
  }
  return false;
}

bool within(Vector3 p, Vector3 r, double rMax, const Grid& g) {
  double r2 = rMax*rMax;
  Vector3 d = g.wrapDiff(p - r);
  return (d.length2() < r2);
}

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if (argc != 10 &&  argc != 11) {
    printf("Make a mask set to 1 for grid nodes within radius of any point.\n");
    printf("Usage: %s srcGrid pointFile paramFile switchDist cutoffDist ewaldFactor Vmax blurCount [x | y | z | xy | xyz] outFile\n", argv[0]);
    return 0;
  }

  const char* srcFile = argv[1];
  const char* pointFile = argv[2];
  const char* paramFile = argv[3];
  const double switchDist = strtod(argv[4], NULL);
  const double cutoffDist = strtod(argv[5], NULL);
  const double ewaldFactor = strtod(argv[6], NULL);
  const double Vmax = strtod(argv[7], NULL);
  const int blurCount = atoi(argv[8]);
  const char* wrap = (argc == 11) ? argv[9] : "";
  const char* outFile = argv[argc-1];
  
  const double convertFactor = 557.0; // convert to kT/e for T = 300K (to match PMEpot plugin)
  
  // Figure out wrapping situation
  bool wrapX = false, wrapY = false, wrapZ = false;
  for (unsigned int i = 0; i < strlen(wrap); i++) {
      wrapX = wrapX || (tolower(wrap[i]) == 'x');
      wrapY = wrapY || (tolower(wrap[i]) == 'y');
      wrapZ = wrapZ || (tolower(wrap[i]) == 'z');
  }
  printf("Wrapping X Y Z = %d %d %d\n", wrapX, wrapY, wrapZ);
  
  // Load the grids.
  Grid src(srcFile);
  src.zero();
  int n = src.length();
  printf("Loaded `%s' which contains %d nodes.\n", srcFile, n); 
  
  // Load the coordinates.
  printf("Loading the coordinates.\n");
  Scatter pos(pointFile);
  printf("Loaded %d points from `%s'.\n", pos.length(), pointFile);
  
  // Load ELEC parameters
  printf("Loading the source parameters.\n");
  Scatter param(paramFile);
  printf("Loaded %d points from `%s'.\n", param.length(), paramFile);
  
  // Generate the cell decomposition.
  double cutoff = 1.2*cutoffDist;
  Matrix3 sys = src.getBox();
  Vector3 origin = src.getOrigin();
  printf("\nGenerating cell decomposition for %s with cutoff = %.10g...\n", argv[1], cutoff);
  CellDecomposition cell(sys, origin, cutoff);
  printf("Cell decomposition of %i cells\n", cell.length());
  cell.decompose(pos);
  printf("Cell decomposition contains %i points\n", cell.countPoints());
  
  printf("cutoffDist = %f\n", cutoffDist);
  
  // Make the mask.
  Grid mask(src);
  int counter = 0;
  int complete = 0;
#pragma omp parallel for
  for (int i = 0; i < n; i++) {
    Vector3 r = src.getPosition(i);
    // Get the indices (in pos) of possible neighboring points.
    IndexList neigh = cell.neighborhood(r);
    const int nNeighs = neigh.length();

    for (int j = 0; j < nNeighs; j++) {
      // Get the position of each possible neighboring point.
      Vector3 p = pos.get(neigh.get(j));
      if (within(p, r, cutoffDist, src)) {
	// wrap selected directions
	Vector3 d1 = p - r;
 	Vector3 d2 = src.wrapDiff(p - r);
	Vector3 d = Vector3(wrapX ? d2.x : d1.x,
			    wrapY ? d2.y : d1.y,
			    wrapZ ? d2.z : d1.z);
	
	double R = d.length();
	Vector3 par = param.get(neigh.get(j));
	double charge = par.x;
	double V = mask.getValue(i);
	
	double V_curr = convertFactor * erfc(ewaldFactor * R) * charge / R; // screened charge
	
	if (R >= switchDist) {
	    double RR = (R - switchDist)/(cutoffDist - switchDist);
	    V_curr *= 1 + RR*RR * (2*RR - 3);
	}
	
	double V_new = (V + V_curr > Vmax) ? Vmax : (V + V_curr < -Vmax) ? -Vmax : V + V_curr;
	//double V_new = V + V_curr;
	mask.setValue(i, V_new);
	//break;
      }
    }

    // Write the progress.
#pragma omp critical
    {
	int comp = (100*(long)counter)/n;
	if (abs(complete - comp) >= 5) {
	    printf("%d percent complete\n", comp);
	    complete = comp;
	}
	counter++;
    }
  }

  // Blur.
  for (int j = 0; j < blurCount; j++) mask.blur();

  char comments[256];
  sprintf(comments, "%s patched", srcFile);
  mask.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
