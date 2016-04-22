///////////////////////////////////////////////////////////////////////  
// Author: Jeff Comer <jcomer2@illinois.edu>    
#include <stdlib.h>
#include <time.h>

#include "useful.H"
#include "Grid.H"
#include "BrownTownGated.H"
#include "Random.H"

String makePdbLine(const String& tempLine, int index, const char* segName, int resId, Vector3 r, double beta) {
  char s[128];

  String record("ATOM  ");
  sprintf(s, "     %5i ", index);
  String si = String(s).range(-6,-1);
  String temp0 = tempLine.range(12,21);
  
  sprintf(s, "    %d", resId);
  String res = String(s).range(-4,-1);
  String temp1 = tempLine.range(26,29);
  
  sprintf(s,"       %.3f", r.x);
  String sx = String(s).range(-8,-1);
  sprintf(s,"       %.3f", r.y);
  String sy = String(s).range(-8,-1);
  sprintf(s,"       %.3f", r.z);
  String sz = String(s).range(-8,-1);

  String temp2 = tempLine.range(54,59);
  sprintf(s,"    %.2f", beta);
  String bet = String(s).range(-6,-1);
  String temp3 = tempLine.range(66,71);

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
  ret.add(bet);
  ret.add(temp3);
  ret.add(seg);
  
  return ret;
}

void writePdbTraj(const char* fileName, const Vector3* atomPos, int atomNum, Vector3 sysDim, double beta) {
  char s[128];
  sprintf(s, "CRYST1   %.3f   %.3f   %.3f  90.00  90.00  90.00 P 1           1\n", sysDim.x, sysDim.y, sysDim.z);
  
  String sysLine(s);
  String tempLine("ATOM      1  CA  MOL S   1      -6.210  -9.711   3.288  0.00  0.00      ION");
  String line;

  FILE* out = fopen(fileName, "w");
  fprintf(out, sysLine.val());

  for (int i = 0; i < atomNum; i++) {
    line = makePdbLine(tempLine, i, "ION", i, atomPos[i], beta);
    fprintf(out, line.val());
    fprintf(out, "\n");
  }
  fprintf(out, "END\n");
  fclose(out);
}

void appendPdbTraj(const char* fileName, const Vector3* atomPos, int atomNum, Vector3 sysDim, double beta) {
  String tempLine("ATOM      1  CA  MOL S   1      -6.210  -9.711   3.288  0.00  0.00      ION");
  String line;

  FILE* out = fopen(fileName, "a");
  for (int i = 0; i < atomNum; i++) {
    line = makePdbLine(tempLine, i, "ION", i, atomPos[i], beta);
    fprintf(out, line.val());
    fprintf(out, "\n");
  }
  fprintf(out, "END\n");
  fclose(out);
}

int main(int argc, char* argv[])
{
  if ( argc != 11 ) {
    printf("Usage: %s gridFile gateGridFile steps diffusion dt externalForceZ num outPeriod gateZ outName\n", argv[0]);
    printf("You entered %i arguments.\n", argc);
    return 0;
  }
  const double kT = 1.0;

  const char* gridFile = argv[1];
  const char* gateGridFile = argv[2];
  const char* outName = argv[argc-1];

  const long int steps = atol(argv[3]);
  const double diffusion = strtod(argv[4],NULL);
  const double dt = strtod(argv[5],NULL);
  const double externalForceZ = strtod(argv[6],NULL);
  const int num = atoi(argv[7]);
  const int outPeriod = atoi(argv[8]);
  const double gateZ = strtod(argv[9],NULL);

  printf("Brownian Dynamics initiated with command:\n");
  for (int i = 0; i < argc; i++) printf("%s ", argv[i]);
  printf("\n");
  printf("Running %d Brownian Dynamics steps with grid %s, diffusion constant %g, and timestep %g.\n", steps, argv[1], diffusion, dt);
  printf("The external force along z is %f.\n", externalForceZ);
  printf("There are %d particles.\n", num);
  printf("The output period is %d.\n", outPeriod);
  printf("A typical displacement in one output period is %g.\n", sqrt(2*diffusion*dt*outPeriod));
  printf("Periodic boundaries in x, y, and z.\n");
  printf("\nParameter list: \n");
  printf("gridFile %s\n", gridFile);
  printf("gateGridFile %s\n", gateGridFile);
  printf("steps %d\n", steps);
  printf("diffusion %.10g\n", diffusion);
  printf("dt %.10g\n", dt);
  printf("externalForceZ %.10g\n", externalForceZ);
  printf("num %d\n", num);
  printf("outPeriod %d\n", outPeriod);
  printf("gateZ %.10g\n", gateZ);
  printf("outName %s\n\n", outName);

  char outTrajFile[256];
  sprintf(outTrajFile, "traj_%s.pdb", outName);
 
  // Load the potential grid.
  Grid pot(gridFile);
  Grid gate(gateGridFile);
  Vector3 extForce = Vector3(0, 0, externalForceZ);

  // Instantiate the Brownian Dynamics object.
  BrownTown brown(pot, gate, diffusion, kT, dt);
  brown.setPeriodic(1,1,1);
  
  // Seed the random number generator.
  int seed = (unsigned int)time((time_t *)NULL);
  printf("Random number generator seed: %i\n", seed);
  Random randoGen(seed);

  // Set initial conditions.
  Vector3* r = new Vector3[num];
  for (int i = 0; i < num; i++) {
    //r[i] = 0.5*(pot.getOrigin() + pot.getCenter());
    r[i] = pot.getCenter();
  }
  bool gated0 = false;

  Vector3 sysDim = pot.getExtent();
  writePdbTraj(outTrajFile, r, num, sysDim, (gated0)?1.0:0.0); 

  // Run the Brownian Dynamics steps.
  for (long int s = 0; s < steps; s++) {
    double t = dt*s;
    
    int currEvents = 0;
    bool gated = false;
    for (int i = 0; i < num; i++) {
      Vector3 rando = randoGen.gaussian_vector();
      r[i] = brown.step(r[i], extForce, rando, gated0);
      if (fabs(r[i].z) < gateZ) gated = true;
    }
    gated0 = gated;

     if (s % outPeriod == 0) {
	double percent = (100.0*s)/steps;
	printf("%.2f percent complete, step %d\n", percent, s);
      
	appendPdbTraj(outTrajFile, r, num, sysDim, (gated0)?1.0:0.0);    
      }
  }

  delete[] r;

  return 0;
}
