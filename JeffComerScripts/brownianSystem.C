///////////////////////////////////////////////////////////////////////  
// Author: Jeff Comer <jcomer2@illinois.edu>    
#include <stdlib.h>
#include <time.h>

#include "useful.H"
#include "Grid1.H"
#include "BrownTown.H"
#include "Random.h"

String makePdbLine(const String& tempLine, int index, const char* segName, int resId, Vector3 r) {
  char s[128];

  String record("ATOM  ");
  sprintf(s, "     %5i ", index);
  String si = String(s).range(-6,-1);
  String temp0 = tempLine.range(12,21);
  
  sprintf(s, "    %d", resId);
  String res = String(s).range(-4,-1);
  String temp1 = tempLine.range(26,29);
  
  sprintf(s,"       %8.3f", r.x);
  String sx = String(s).range(-8,-1);
  sprintf(s,"       %8.3f", r.y);
  String sy = String(s).range(-8,-1);
  sprintf(s,"       %8.3f", r.z);
  String sz = String(s).range(-8,-1);

  String temp2 = tempLine.range(54,71);
  sprintf(s, "%s    ", segName);
  String seg = String(s).range(0,3);

  String ret(record);
  ret.add(si);
  ret.add(temp0);
  ret.add(res);
  ret.add(temp1);
  ret.add(sx);
  ret.add(sy);
  ret.add(sz);
  ret.add(temp2);
  ret.add(seg);
  
  return ret;
}

void writePdbTraj(const char* fileName, const Vector3* atomPos, int atomNum, Vector3 sysDim) {
  char s[128];
  sprintf(s, "CRYST1   %.3f   %.3f   %.3f  90.00  90.00  90.00 P 1           1\n", sysDim.x, sysDim.y, sysDim.z);
  
  String sysLine(s);
  String tempLine("ATOM      1  CA  MOL S   1      -6.210  -9.711   3.288  0.00  0.00      SIO");
  String line;

  FILE* out = fopen(fileName, "w");
  fprintf(out, sysLine.val());

  for (int i = 0; i < atomNum; i++) {
    line = makePdbLine(tempLine, i, "DMMP", i, atomPos[i]);
    fprintf(out, line.val());
    fprintf(out, "\n");
  }
  fprintf(out, "END\n");
  fclose(out);
}

void appendPdbTraj(const char* fileName, const Vector3* atomPos, int atomNum, Vector3 sysDim) {
  String tempLine("ATOM      1  CA  MOL S   1      -6.210  -9.711   3.288  0.00  0.00      SIO");
  String line;

  FILE* out = fopen(fileName, "a");
  for (int i = 0; i < atomNum; i++) {
    line = makePdbLine(tempLine, i, "DMMP", i, atomPos[i]);
    fprintf(out, line.val());
    fprintf(out, "\n");
  }
  fprintf(out, "END\n");
  fclose(out);
}

int main(int argc, char* argv[])
{
  if ( argc != 9 ) {
    printf("Usage: %s gridFile diffusion dt steps surfZ edgeZ num outName\n", argv[0]);
    printf("%i\n", argc);
    return 0;
  }

  double kT = 1.0;
  double diffusion = strtod(argv[2],NULL);
  double dt = strtod(argv[3],NULL);
  long int steps = atol(argv[4]);
  double surfZ = strtod(argv[5],NULL);
  double edgeZ = strtod(argv[6],NULL);
  int num = atoi(argv[7]);
  printf("Running %d Brownian Dynamics steps with grid %s, diffusion constant %g, and timestep %g.\n", steps, argv[1], diffusion, dt);
  printf("Periodic boundaries in x, y.\n");
  printf("Molecules with z < %g are considered to be bound.\n", surfZ);
  printf("Forces for z > %g assumed to be zero.\n", edgeZ);  
  const char* outName = argv[argc-1];

  int displayPeriod = 1000;
  char outForceFile[256];
  char outEventFile[256];
  char outTrajFile[256];
  char outBoundFile[256];
  sprintf(outForceFile, "fz_%s.dat", outName);
  sprintf(outEventFile, "event_%s.dat", outName);
  sprintf(outTrajFile, "traj_%s.pdb", outName);
  sprintf(outBoundFile, "bound_%s.dat", outName);

  // Load the potential grid.
  Grid pot(argv[1]);

  // Instantiate the Brownian Dynamics object.
  BrownTown brown(pot, diffusion, kT, dt, edgeZ);
  brown.setPeriodic(1,1,0);
  
  // Seed the random number generator.
  int seed = (unsigned int)time((time_t *)NULL);
  printf("Random number generator seed: %i\n", seed);
  Random randoGen(seed);

  // Open the file.s
  FILE* outBound = fopen(outBoundFile, "w");
  FILE* outEvent = fopen(outEventFile, "w");

  // Set initial conditions.
  Vector3* r = new Vector3[num];
  double* eventTime = new double[num];
  bool* event = new bool[num];
  for (int i = 0; i < num; i++) {
    //r[i] = pot.getCenter();
    r[i] = pot.getDestination();
    eventTime[i] = 0.0;
  }
  double boundTime = 0.0;
  Vector3 sysDim = pot.getExtent();
  writePdbTraj(outTrajFile, r, num, sysDim); 

  // Run the Brownian Dynamics steps.
  for (long int s = 0; s < steps; s++) {
    double t = dt*s;
    
    int currEvents = 0;
    for (int i = 0; i < num; i++) {
      Vector3 rando = randoGen.gaussian_vector();
      r[i] = brown.step(r[i], rando);

      // Track binding events.
      if (event[i]) {
	currEvents++;
	if (r[i].z > surfZ) {
	  // A complete binding event.
	  double dur = t - eventTime[i];
	  fprintf(outEvent, "%.10g %.10g\n", eventTime[i], dur);
	  event[i] = false;
	  boundTime += dur;
	}
      } else {
	if (r[i].z < surfZ) {
	  // New binding event
	  eventTime[i] = t;
	  event[i] = true;
	}
      }
    }
    
    
    if (s % displayPeriod == 0) {
      int percent = (100*s)/steps;
      printf("%d percent complete, step %d\n", percent, s);
      
      appendPdbTraj(outTrajFile, r, num, sysDim);    
      double currBoundFrac = (double)currEvents/(double)num;
      fprintf(outBound, "%.10g %.10g\n", t, currBoundFrac);
    }
  }
  fclose(outBound);
  fclose(outEvent);

  double boundFrac = boundTime/(dt*steps*num);
  printf("Bound fraction: %.10g\n", boundFrac);
  return 0;
}
