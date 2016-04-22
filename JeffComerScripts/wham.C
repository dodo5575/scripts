///////////////////////////////////////////////////////////////////////
// Use the weighted histogram analysis method.
// Author: Jeff Comer <jcomer2@illinois.edu>

#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "useful.H"
#include "Grid.H"
#include "Scatter.H"
#include "HistogramGrid.H"

class UmbrellaWindow {
public:
  UmbrellaWindow() {
    steps = 0;
  }

  UmbrellaWindow(const Scatter& source, Vector3 center0, Vector3 spring0) {
    steps = source.length();
    center = center0;
    spring = spring0;
  }

public:
  Vector3 center;
  Vector3 spring;
  int steps;
};

class BinDim {
public:
  BinDim(Matrix3 box, Vector3 origin0, int nx0, int ny0, int nz0) {
    basis = Matrix3(box.ex()/nx0, box.ey()/ny0, box.ez()/nz0);
    origin = origin0;
    nx = nx0;
    ny = ny0;
    nz = nz0;
  }
  
  BinDim(Vector3 r0, Vector3 r1, int nx0, int ny0, int nz0) {
    origin = r0;
    Vector3 d = r1 - r0;
    basis = Matrix3(d.x/nx0, d.y/ny0, d.z/nz0);
    nx = nx0;
    ny = ny0;
    nz = nz0;
  }

  BinDim(const Scatter** win, int winN, int nx0, int ny0, int nz0) {
    Vector3 b;
    Vector3 boundMin = win[0]->minBound();
    Vector3 boundMax = win[0]->maxBound();
    for (int w = 1; w < winN; w++) {
      b = win[w]->minBound();
      if (b.x < boundMin.x) boundMin.x = b.x;
      if (b.y < boundMin.y) boundMin.y = b.y;
      if (b.z < boundMin.z) boundMin.z = b.z;

      b = win[w]->maxBound();
      if (b.x > boundMax.x) boundMax.x = b.x;
      if (b.y > boundMax.y) boundMax.y = b.y;
      if (b.z > boundMax.z) boundMax.z = b.z;
    }

    // Don't allow bins of zero size.
    if (boundMax.x - boundMin.x == 0.0) {
      printf("Expanding very small bin along x.\n");
      boundMax.x += 1.0;
      boundMin.x -= 1.0;
    } 
    if (boundMax.y - boundMin.y == 0.0) {
      printf("Expanding very small bin along y.\n");
      boundMax.y += 1.0;
      boundMin.y -= 1.0;
    }
    if (boundMax.z - boundMin.z == 0.0) {
      printf("Expanding very small bin along z.\n");
      boundMax.z += 1.0;
      boundMin.z -= 1.0;
    }

    nx = nx0;
    ny = ny0;
    nz = nz0;
    origin = boundMin;
    Matrix3 box = Matrix3(boundMax.x-boundMin.x, boundMax.y-boundMin.y, boundMax.z-boundMin.z); 
    basis = Matrix3(box.ex()/nx0, box.ey()/ny0, box.ez()/nz0);
  }

  int length() const {
    return nx*ny*nz;
  }

public:
  int nx, ny, nz;
  Vector3 origin;
  Matrix3 basis;
};

double springEnergy(Vector3 spring, Vector3 r, Vector3 r0) {
  Vector3 d = r - r0;
  return 0.5*(spring.x*d.x*d.x + spring.y*d.y*d.y + spring.z*d.z*d.z);
}

void writeArray(const char* outFile, double* d, int n) {
  FILE* out = fopen(outFile, "w");
  for (int i = 0; i < n; i++) {
    fprintf(out, "%.10g\n", d[i]);
  } 
}

// Compute the PMF from a probability distribution.
// Write the result and a mask with the valid points.
void writePmf(const char* fileName, const char* comments, const Grid& prob) {
  int b;
  const int n = prob.getSize();

  char maskName[256];
  sprintf(maskName, "valid_%s", fileName);

  // Find the minimum and maximum of the PMF.
  double p = 0.0;
  for (b = 0; b < n; b++) {
    p = prob.getValue(b);
    if (p > 0.0) break;
  }
  if (b == n) {
    printf("Error! Probability distribution is zero.\n");
    return;
  }

  double pmfMin = -log(p);
  double pmfMax = -log(p);
  for (b = 0; b < n; b++) {
    if (prob.getValue(b) > 0.0) {
      double u = -log(prob.getValue(b));
      if (u < pmfMin) pmfMin = u;
      if (u > pmfMax) pmfMax = u;
    }
  }
  
  // Calculate the PMF.
  Grid pmf(prob);
  Grid mask(prob);
  // The PMF is associated with the bin centers.
  //pmf.shiftToCenters(); 
  //mask.shiftToCenters();
  for (b = 0; b < n; b++) {
    if (prob.getValue(b) > 0.0) {
      double u = -log(prob.getValue(b));
      pmf.setValue(b, u - pmfMin);
      mask.setValue(b, 1);
    } else {
      pmf.setValue(b, pmfMax - pmfMin);
      mask.setValue(b, 0);
    }
  }

  pmf.write(fileName, comments);
  mask.write(maskName, comments);

  char profName[256];
  sprintf(profName, "%s.x.dat", fileName);
  pmf.averageProfile(profName, 0);
  sprintf(profName, "%s.y.dat", fileName);
  pmf.averageProfile(profName, 1);
  sprintf(profName, "%s.z.dat", fileName);
  pmf.averageProfile(profName, 2);
}

// Omit windows with no data or large average energies.
// Deallocate their memory.
int filterWindows(Scatter** sim, UmbrellaWindow* win, int winN0, double maxWinEnergy) {
  int wg = 0;
  for (int w = 0; w < winN0; w++) {
    // Compute the average energy for each window.
    double energy = 0.0;
    Vector3 pos(0.0);
    for (int i = 0; i < sim[w]->length(); i++) {
      energy += springEnergy(win[w].spring, win[w].center, sim[w]->get(i));
      pos += sim[w]->get(i);
    }
    energy /= sim[w]->length();
    pos /= sim[w]->length();

    // Omit windows that have an energy that is too large.
    if (energy < maxWinEnergy && sim[w]->length() > 0) {
      sim[wg] = sim[w];
      win[wg] = win[w];
      wg++;
      printf("WINDOW: %i %s ", w, win[w].center.toString().val());
    }  else {
      printf("OMIT WINDOW: %i %s ", w, win[w].center.toString().val());
      delete sim[w];
    }
    printf("MEAN_ENERGY: %g MEAN_POSITION: %s\n", energy, pos.toString().val());
  }
  const int winN = wg;

  printf("%d of %d windows contribute.\n", winN, winN0);
  return winN;
}

// Make histogram for each window and put it in hist. Windows that do not contribute to the histogram will be omitted.
void histogramWindows(HistogramGrid** hist, const Scatter** sim, const UmbrellaWindow* win, int winN, int binNx, int binNy, int binNz) {
  int w;
  // Form the bins to cover all included windows.
  BinDim bin(sim, winN, binNx, binNy, binNz);
  printf("bin: %s\n", bin.basis.toString().val());


  // Histogram the data to form the biased probability distributions for each window.
  printf("Histogramming data...\n");
  for (w = 0; w < winN; w++) {
    hist[w] = new HistogramGrid(bin.basis, bin.origin, bin.nx, bin.ny, bin.nz);

    int histCount = hist[w]->add(*sim[w]);
    printf("Histogrammed %d of %d points for window %s.\n", histCount, sim[w]->length(), win[w].center.toString().val());
  }
  printf("Done.\n");
}

// Run the Weighted Histogram Analysis Method in three dimensions.
// winIndex contains the indices of the valid windows in "win".
void wham3d(const char* outName, HistogramGrid** hist, const UmbrellaWindow* win, int winN, double tol) {
  int w, b;

  char outFile[256];
  sprintf(outFile, "pmf_%s.dx", outName);
  char tempFile[256];
  sprintf(tempFile, "preview_%s.dx", outName);
  char histFile[256];
  sprintf(histFile, "hist_%s.dx", outName);
  
  printf("\nComputing the PMF distribution in three dimensions by the WHAM method with %d windows.\n", winN);
  printf("Based on B. Roux. Computer Physics Communications 91, 275 (1995).\n");
  
  // Get the histogram counts.
  double* histCount = new double[winN];
  for (w = 0; w < winN; w++) histCount[w] = hist[w]->total();

  // Write the total histogram.
  Grid total(*hist[0]);
  total.zero();
  for (w = 0; w < winN; w++) total.addGrid(*hist[w]);
  //total.shiftToCenters();
  total.write(histFile, "total histogram");
  printf("Wrote the total histogram to %s.\n", histFile);
 
  // Scale the histograms to make them probability distributions.
  const int binN = hist[0]->length();
  const double binVol = hist[0]->getCellVolume();  
  for (w = 0; w < winN; w++) hist[w]->scale(1.0/(histCount[w]*binVol));
  
  // Compute the numerator of Roux Eq. 8 beforehand.
  printf("\nPrecomputing the numerator of Roux Eq. 8...\n");
  double* binNumer = new double[binN];
  for (b = 0; b < binN; b++) {
    double numer = 0.0;
    
    for(w = 0; w < winN; w++) {
      // Get the biased prob. distrib. at this bin center for this window.
      double p = hist[w]->getValue(b);

      // Add to the numerator sum in Roux Eq. 8.
      numer += histCount[w]*p;
    }
    binNumer[b] = numer;
  }
  printf("Done.\n");
  //writeArray("numer_C.txt", binNumer, binN);

  // Compute exp(-uBias) beforehand for Roux Eq. 8 and 9.
  printf("\nPrecomputing exp(-uBias) for Roux Eq. 8 and 9...\n");
  double** biasFactor = new double*[binN];
  for (b = 0; b < binN; b++) {
    biasFactor[b] = new double[winN];
    Vector3 binCen = hist[0]->getBinCenter(b);
    for (w = 0; w < winN; w++) {
      double uBias = springEnergy(win[w].spring, win[w].center, binCen);
      biasFactor[b][w] = exp(-uBias);
      //double d = (binCen-win[w].center).length();
    }
  }
  printf("Done.\n");

  // Make the initial guess for the f constants.
  double* winF = new double[winN];
  double* winFOld = new double[winN];
  for (w = 0; w < winN; w++) {
    winF[w] = 0.0;
    winFOld[w] = 0.0;
  }

  // This contains the estimate of the unbiased prob. distrib.
  Grid prob(*hist[0]);
  prob.zero();

  //////////////////////////////
  //
  // Iterate using the WHAM equations.
  printf("\nBeginning WHAM iterations.\n");
  double maxChange = 1.0;
  int iter = 0;
  while (maxChange > tol) {
    //////////////////////////////
    // Roux Eq. 8.
    //
    // Estimate the unbiased probability distribution.
    for (b = 0; b < binN; b++) {
      double numer = binNumer[b];
      double denom = 0.0;
      // Loop through the windows to do the sum.
      for (w = 0; w < winN; w++)
	denom += histCount[w]*exp(winF[w])*biasFactor[b][w];
      prob.setValue(b, numer/denom);
    }
    
    //////////////////////////////
    // Roux Eq. 9
    //
    // Integrate to obtain the new f constants for each window.
    double* ptr = winFOld;
    winFOld = winF;
    winF = ptr;
    double sum;
    for (w = 0; w < winN; w++) {
      // Numerically integrate over the reaction coordinates (over the bins).
      sum = 0.0;
      for (b = 0; b < binN; b++)
	sum += binVol*biasFactor[b][w]*prob.getValue(b);
      winF[w] = -log(sum);
    }
 
    // Check for convergence.
    w = 1;
    double dfOld = winFOld[w] - winFOld[w-1];
    double df = winF[w] - winF[w-1];
    maxChange = fabs(df-dfOld);
    for (w = 2; w < winN; w++) {
      dfOld = winFOld[w] - winFOld[w-1];
      df = winF[w] - winF[w-1];
      if (fabs(df-dfOld) > maxChange) {
	maxChange = fabs(df-dfOld);
      }
    }
    printf("WHAM iteration %d: %.5g\n", iter, maxChange);
    iter++;

    // Write the intermediate results.
    if (iter % 100 == 0) {
      char comm[256];
      sprintf(comm, "PMF iteration %d", iter);
      writePmf(tempFile, comm, prob);
    }
  }

  // Write the F constants.
  char outNameF[256];
  sprintf(outNameF, "fconst_%s.dat", outName);
  FILE* outF = fopen(outNameF, "w");
  for (w = 0; w < winN; w++) fprintf(outF, "%.10g\n", winF[w]);
  fclose(outF);

  // We are finished with the WHAM iterations.
  printf("\nWHAM iterations finished.\n");

  // Write the PMF.
  writePmf(outFile, "PMF", prob);
  printf("Wrote %s.\n", outFile);

  // Clean up the memory.
  for (b = 0; b < binN; b++)  delete biasFactor[b];

  delete[] binNumer;
  delete[] biasFactor;
  delete[] winF;
  delete[] winFOld;
}

void readIndexFile(const char* fileName, int num, String* prefix, UmbrellaWindow* win) {
  int nRead;
  int n = 0;
  double x, y, z, kx, ky ,kz;
  char line[256];
  char ind[256];

  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("readIndexFile Couldn't open file %s\n.",fileName);
    exit(-1);
  }

  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead=sscanf(line,"%s %lf %lf %lf %lf %lf %lf",ind,&x,&y,&z,&kx,&ky,&kz);
    if (nRead >= 7) {
      prefix[n].add(ind);
      win[n].center.x = x;
      win[n].center.y = y;
      win[n].center.z = z;
      win[n].spring.x = kx;
      win[n].spring.y = ky;
      win[n].spring.z = kz;
      n++;
      if (n >= num) break;
    } else {
      printf("Warning! Improperly formatted index file %s.\n", fileName);
      printf("Simulation index file format: winDataFile x0 y0 z0 kx ky kz.\n");
    }
  }
    
  fclose(inp);
}

int countIndexFile(const char* fileName) {
  int nRead;
  int n = 0;
  double x, y, z, kx, ky, kz;
  char line[256];
  char ind[256];

  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("countIndexFile Couldn't open file %s\n.",fileName);
    exit(-1);
  }

  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead=sscanf(line,"%s %lf %lf %lf %lf %lf %lf",ind,&x,&y,&z,&kx,&ky,&kz);
    if (nRead >= 7) n++;
  }
    
  fclose(inp);
  return n;
}

void readBinDimensions(const char* fileName, int* pnx, int* pny, int* pnz) {
  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("readBinDimensions Couldn't open file %s\n.",fileName);
    exit(-1);
  }
  int nRead;

  char line[256];
  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead=sscanf(line, "%d %d %d",pnx,pny,pnz);
    if (nRead != 3) {
      printf("readBinDimensions File %s improperly formatted.\n", fileName);
      printf("Format is 'nx ny nz'.\n");
      exit(1);
    }
  }

  return;
}

///////////////////////////////////////////////////////////////////////
// Drivers
// indexFile has records "prefix x y z kx ky kz".
// Spring constants are converted using the given kT.
// Output is in kT.
int main(int argc, char* argv[]) {
  if ( argc != 6 ) {
    printf("Usage: %s indexFile (binDimensionsFile | gridFile.dx) cutTime0 kT outName\n", argv[0]);
    return 0;
  }

  // Parameters:
  double tol = 0.0001;
  //Vector3 springK(0.8, 0.8, 4.0); // in kcal/mol A^-2
  const double cutTime0 = strtod(argv[3], NULL);
  const double kBT = strtod(argv[4], NULL);
  // Input:
  const char* indexFile = argv[1];
  const char* binDimFile = argv[2];
  // Output:
  const char* outName = argv[5];
  
  //const double kBT = 0.5862292; // in kcal/mol at 295 K
  const double maxWinEnergy = 20.0; // Windows with average energies greater than this are discarded.
  int w;

  
  

  int winN = countIndexFile(indexFile);
  printf("Found %d windows in simulation index file %s.\n", winN, indexFile);
  
  // Window variables
  String* winDataFile = new String[winN];
  UmbrellaWindow* win = new UmbrellaWindow[winN];
  readIndexFile(indexFile, winN, winDataFile, win);

  // Convert to k_B T/A^2.
  for (w = 0; w < winN; w++) win[w].spring /= kBT;
  printf("Converted the spring constants from E^2/x^2 to kT/x^2 using E = %g kT.\n", kBT);

  // Read the simulation data.
  printf("Reading simulation data, ignoring data with time < %g.\n", cutTime0);
  Scatter** sim = new Scatter*[winN];
  for (w = 0; w < winN; w++) {
    sim[w] = new Scatter(winDataFile[w], cutTime0);
    win[w].steps = sim[w]->length();

    printf("Read %d data points from %s.\n", win[w].steps, winDataFile[w].val());
  }

  // Filter the windows.
  const int winN1 = filterWindows(sim, win, winN, maxWinEnergy);

  // Read the bin dimensions from a file.
  Grid* binDimGrid;
  bool isPeriodic; 
  int len = strlen(binDimFile);
  if (strcmp(binDimFile + len-3, ".dx") == 0) {
    printf("The binDimensions file is a dx file (`%s'). The data will be wrapped.\n", binDimFile);
    
    binDimGrid = new Grid(binDimFile);
    isPeriodic = true;
  } else {
    printf("Read bin dimensions from %s.\n", binDimFile);

    binDimGrid = new Grid()
    isPeriodic = false;
  }

  int binNx, binNy, binNz; 
  readBinDimensions(binDimFile, &binNx, &binNy, &binNz);

  // Make the histograms.
  HistogramGrid** hist = new HistogramGrid*[winN1];
  histogramWindows(hist, (const Scatter**)sim, win, winN1, binNx, binNy, binNz);

  

  // Run WHAM.
  wham3d(outName, hist, win, winN1, tol);

  printf("The spring constants were converted to kT/length^2 using the factor %g.\n", kBT);
  printf("The output is in kT.\n");

  // Clean up.
  for (w = 0; w < winN1; w++) delete sim[w];
  for (w = 0; w < winN1; w++)  delete hist[w];
 
  delete[] sim;
  delete[] win;
  delete[] winDataFile;
  delete[] hist;

  return 0;
}
