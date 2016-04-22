///////////////////////////////////////////////////////////////////////  
// Author: Jeff Comer <jcomer2@illinois.edu>
#include <ctime>
#include <cstdio>
#include <cstdlib>

#include "useful.H"
#include "Grid.H"
#include "Reader.H"
#include "RandomGsl.H"
#include "ComputeForce.H"
#include "BrownianParticleType.H"

using namespace std;

class GrandBrownTown {
public:
  GrandBrownTown(const char* configFile, const char* outArg, const long int randomSeed) {
    sprintf(outTrajFile, "%s.traj", outArg);
    sprintf(outCurrFile, "curr_%s.dat", outArg);
    sprintf(outPdbFile, "%s.pdb", outArg);

    // Read the parameters.
    setDefaults();
    readParameters(configFile);

    // Load the potential grids.
    printf("Loading the potential grids...\n");
    for (int i = 0; i < numParts; i++) {
      part[i].grid = new Grid(partGridFile[i].val());
      printf("Loaded %s.\n", partGridFile[i].val());
      printf("System size %s.\n", part[i].grid->getExtent().toString().val());
    }
    // Get the system dimensions.
    sysDim = part[0].grid->getExtent();
    
    // Seed the random number generator. 
    long int randoSeed;
    if (seed == 0) {
      long int randomSeed1 = randomSeed;
      for (int i = 0; i < 4; i++) randomSeed1 *= randomSeed1 + 1;
      randoSeed = (unsigned int)time((time_t *)NULL) + randomSeed1;
    } else randoSeed = seed + randomSeed;
    randoGen.init(randoSeed);
    printf("\nRandom number generator seed: %ld\n", randoSeed);

    // Get the total number of particles.
    num0 = 0; // mean number of particles
    for (int i = 0; i < numParts; i++) num0 += part[i].num;
    numCap = numCapFactor*num0; // max number of particles
    num = num0;
    printf("Initial number of particles: %d\n", num);

    // Allocate particle variables.
    pos = new Vector3[numCap];
    pos1 = new Vector3[numCap];
    type = new int[numCap];
    posLast = new Vector3[numCap];

    // Set the initial conditions.
    printf("\n");
    Vector3 buffer = (part[0].grid->getCenter() + 2.0*part[0].grid->getOrigin())/3.0;
    initialZ = buffer.z;
    populate();
    if (inputCoordinates.length() == 0) {
      initialCond(); 
      printf("Set random initial conditions.\n");
    } else {
      bool loaded = loadCoordinates(inputCoordinates.val());
      if (loaded) printf("Loaded initial coordinates from %s.\n", inputCoordinates.val());
      else {
	printf("WARNING: Could not read input coordinate file %s.\n",  inputCoordinates.val());
	initialCond(); 
	printf("Set random initial conditions.\n");
      }
    }

    // Prepare for internal force computation.
    printf("\n");
    switchStart = cutoff - switchLen;
    if (fullElect == 0)
      printf("Cutting off the potential from %.10g to %.10g.\n", switchStart, switchStart+switchLen);
    internal = new ComputeForce(num, part, numParts, *(part[0].grid), switchStart, switchLen, coulombConst);
    forceInternal = new Vector3[numCap];
    for (int i = 0; i < num; i++) forceInternal[i] = 0.0;
    printf("Prepared for internal force computation.\n");
    internal->decompose(pos);
    IndexList decompDim = internal->decompDim();
    printf("Created a cell decomposition with cutoff of %.10g and dimensions %s.\n", internal->decompCutoff(), decompDim.toString().val());

    // Load the tabulated forces.
    if (tabulatedPotential) {
      printf("Loading the tabulated potentials...\n");
      for (int p = 0; p < numParts*numParts; p++)
	if (partTableFile[p].length() > 0) {
	  int type0 = partTableIndex0[p];
	  int type1 = partTableIndex1[p];

	  internal->addTabulatedPotential(partTableFile[p].val(), type0, type1);
	  printf("Loaded %s for types %s and %s.\n", partTableFile[p].val(), 
		 part[type0].name.val(), part[type1].name.val());
	}
    }
  }

  ~GrandBrownTown() {
    delete[] pos;
    delete[] pos1;
    delete[] posLast;
    delete[] type;
    delete[] part;
    delete[] partGridFile;
    delete[] partTableFile;
    delete[] partTableIndex0;
    delete[] partTableIndex1;
    delete[] forceInternal;
    delete internal;
  }

  Vector3 step(Vector3 r0, Vector3 force, double diffusion) {
    // Get the random kick.
    Vector3 rando = randoGen.gaussian_vector();

    // Perform the Brownian Dynamics step.
     Vector3 r = r0 + force*diffusion*timestep/kT + sqrt(2*diffusion*timestep)*rando;

     // Apply the boundary conditions.
     return wrap(r);
  }

  ////////////////////////////////////////////////////////////////
  // Run the Brownian Dynamics steps.
  void run() {
    printf("\n");
    clock_t clockInit = clock();
    clock_t clock0 = clockInit;
    clock_t clock1;
    newCurrent();
    if (outputPdb == 1) newPdbTraj(outPdbFile);
    else newTrajectory(0.0);

    remember(0.0);
    for (long int s = 1; s < steps; s++) {
      // Compute the internal forces.
      double energy = 0.0;
      if (interparticleForce) {
	if (tabulatedPotential) {
	  switch (fullElect) {
	  case 0:
	    // Remake the cell decomposition.
	    if (s % decompPeriod == 0) internal->decompose(pos);
	    // Compute using the cell decomposition.
	    energy = internal->computeTabulated(forceInternal, pos, type);
	    break;
	  default:
	    // Compute long range electrostatic forces.
	    energy = internal->computeTabulatedFull(forceInternal, pos, type);
	    break;
	  }
	} else {
	  switch (fullElect) {
	  case 0:
	    // Remake the cell decomposition.
	    if (s % decompPeriod == 0) internal->decompose(pos);
	    // Compute using the cell decomposition.
	    energy = internal->compute(forceInternal, pos, type);
	    break;
	  case 1:
	    // Compute long range electrostatic forces.
	    energy = internal->computeFull(forceInternal, pos, type);
	    break;
	  case 2:
	    // Compute only softcore forces.
	    energy = internal->computeSoftcoreFull(forceInternal, pos, type);
	    break;
	  }
	}
      }
    
      // Loop through the particles.
      for (int i = 0; i < num; i++) {
	// Compute the external forces.
	Vector3 forceExternal = Vector3(0.0, 0.0, part[type[i]].charge*electricField);
	Vector3 forceGrid = part[type[i]].grid->interpolateForce(pos[i]);
	//Vector3 forceGrid = 0.0;

	// Compute the total force.
	Vector3 force = forceInternal[i] + forceExternal + forceGrid;

	// Step
	pos1[i] = step(pos[i], force, part[type[i]].diffusion);
      }

       
      // Swap the position pointers.
      Vector3* temp = pos;
      pos = pos1;
      pos1 = temp;
      
      if (s % outputPeriod == 0) {
	double t = s*timestep;

	// Write the current.
	if (currentSegmentZ <= 0.0) writeCurrent(t);
	else writeCurrentSegment(t, currentSegmentZ);

	if (outputPdb) writePdbTraj(outPdbFile);
	else writeTrajectory(t);

	// Handle particle fluctuations.
	if (numberFluct == 1) numberFluctuation();
	
	// Store the current positions.
	remember(t);
	
	if (s % (10*outputPeriod) == 0) {
	  printf("ENERGY: %.10g %.10g\n", t, energy);
	  
	  double percent = (100.0*s)/steps;
	  clock1 = clock();
	  double stepTime = diffClock(clock0, clock1)/steps;
	  printf("step %ld, time %g, %.2f percent complete, %.10g ms/step, %d particles\n", s, t, percent, stepTime, num);
	  clock0 = clock1;
	}
      }
    }
    
    printf("Run time: %.10g s\n", diffClock(clockInit, clock())/1000.0);
  }

  // Get the $n farthest particles in $ind.
  IndexList farthest(IndexList& ind, int n) {
    double* dist = new double[ind.length()];
    for (int i = 0; i < ind.length(); i++) dist[i] = fabs(pos[ind.get(i)].z);

    int* farth = new int[n];
    for (int i = 0; i < n; i++) farth[i] = i;
    insertionSort(farth, dist, n);

    for (int i = n; i < ind.length(); i++) {
      if(dist[i] > dist[farth[0]]) {
	// Substitute the lowest list item.
	farth[0] = i;
	// Resort the list.
	insertionSort(farth, dist, n);
      }
    }
    //for (int i = 0; i < n; i++) printf("%d %g\n", farth[i], dist[farth[i]]);

    IndexList ret;
    for (int j = 0; j < n; j++) ret.add(ind.get(farth[j]));
    //printf("ind: %s\n", ind.toString().val());
    //printf("farth: %s\n", ret.toString().val());

    delete[] farth;
    delete[] dist;

    return ret;
  }

  void setDefaults() {
    // Global parameters
    outputName = "out";
    timestep = 1e-4;
    steps = 100;
    seed = 0;
    inputCoordinates = "";
    numberFluct = 0;
    numberFluctPeriod = 200;
    interparticleForce = 1;
    tabulatedPotential = 0;
    fullElect = 1;
    kT = 1.0;
    coulombConst = 566.440698/92.0;
    electricField = 0.0;
    cutoff = 10.0;
    outputPeriod = 200;
    outputPdb = 0;
    currentSegmentZ = -1.0;

    // Hidden parameters
    decompPeriod = 5;
    switchLen = 2.0;
    maxInitialPot = 0.5;
    numCapFactor = 10;
  }

  void populate() {
    int pn = 0;
    int p = 0;

    for (int i = 0; i < num; i++) {
      type[i] = p;

      pn++;
      if (pn >= part[p].num) {
	p++;
	pn = 0;
      }
    }
  }

  Vector3 findPos(int typ) {
    Vector3 r;
    do {
      r = wrap(Vector3(sysDim.x*randoGen.uniform(),sysDim.y*randoGen.uniform(), sysDim.z*randoGen.uniform()));
    } while (part[typ].grid->interpolatePotential(r) > maxInitialPot);
    return r;
  } 

  Vector3 findPos(int typ, double minZ) {
    Vector3 r;
    do {
      r = wrap(Vector3(sysDim.x*randoGen.uniform(),sysDim.y*randoGen.uniform(), sysDim.z*randoGen.uniform()));
    } while (part[typ].grid->interpolatePotential(r) > maxInitialPot && fabs(r.z) > minZ);
    return r;
  }
 
  bool loadCoordinates(const char* fileName) {
    char line[256];
    FILE* inp = fopen(fileName, "r");
    
    int count = 0;
    while (fgets(line, 256, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;

      String s(line);
      int numTokens = s.tokenCount();
      if (numTokens != 3) {
	printf("ERROR: Invalid coordinate file line: %s\n", line);
	fclose(inp);	
	return false;
      }

      String* tokenList = new String[numTokens];
      s.tokenize(tokenList);
      if (tokenList == NULL) {
	printf("ERROR: Invalid coordinate file line: %s\n", line);
	fclose(inp);
	return false;
      }
            
      if (count >= num) {
	printf("ERROR: Too many coordinates in coordinate file.\n");
	fclose(inp);
	return false;
      }
      double x = strtod(tokenList[0],NULL);
      double y = strtod(tokenList[1],NULL);
      double z = strtod(tokenList[2],NULL);
      pos[count] = Vector3(x,y,z);
      pos1[count] = pos[count];
      count++;

      delete[] tokenList;
    }
    fclose(inp);

    if (count < num) {
      printf("ERROR: Too few coordinates in coordinate file.\n");
      return false;
    }
    return true;
  }

  void initialCondCen() {
    for (int i = 0; i < num; i++) {
      pos[i] = part[type[i]].grid->getCenter();
      pos1[i] = part[type[i]].grid->getCenter();
    }
  }

  void initialCond() {
    for (int i = 0; i < num; i++) {
      // Get the initial positions.
      pos[i] = findPos(type[i], initialZ);
      pos1[i] = pos[i];
    }
  }

  static void insertionSort(int* item, double* key, int n) {
    int i, j, low;

    for (i = 1; i < n; i++) {
      low = item[i];
      j = i;
      while ((j > 0) && (key[item[j-1]] > key[low])) {
	item[j] = item[j-1];
	j--;
      }
      item[j] = low;
    }
  }

  bool readTableFile(const String& value, int currTab) {
    int numTokens = value.tokenCount('@');
    if (numTokens != 3) {
      printf("ERROR: Invalid tabulatedFile: %s\n", value.val());
      return false;
    }

    String* tokenList = new String[numTokens];
    value.tokenize(tokenList, '@');
    if (tokenList == NULL) {
      printf("ERROR: Invalid tabulatedFile: %s\n", value.val());
      return false;
    }

    partTableIndex0[currTab] = atoi(tokenList[0]);
    partTableIndex1[currTab] = atoi(tokenList[1]);
    partTableFile[currTab] = tokenList[2];

    printf("part: %d %d %s\n", partTableIndex0[currTab], partTableIndex1[currTab], partTableFile[currTab].val() );

    return true;
  }

  static double wrap(double x, double l) {
    int image = int(floor(x/l));
    x -= image*l;
    
    if (x > 0.5*l) x -= l;
    return x;
  }

  Vector3 wrap(Vector3 r) const {
    r.x = wrap(r.x, sysDim.x);
    r.y = wrap(r.y, sysDim.y);
    r.z = wrap(r.z, sysDim.z);
    return r;
  }

  double diffClock(clock_t clock1,clock_t clock2)
  {
    double diffTicks = double(clock2-clock1);
    double diffms = (diffTicks*1000)/CLOCKS_PER_SEC;
    return diffms;
  } 

 void newTrajectory(double t) const {
    FILE* out = fopen(outTrajFile, "w");
    for (int i = 0; i < num; i++)
      fprintf(out, "%d %.10g %.10g %.10g %.10g\n", type[i], t, pos[i].x, pos[i].y, pos[i].z);
    fprintf(out, "END\n");
    fclose(out);
  }

  void writeTrajectory(double t) const {
    FILE* out = fopen(outTrajFile, "a");
    for (int i = 0; i < num; i++)
      fprintf(out, "%d %.10g %.10g %.10g %.10g\n", type[i], t, pos[i].x, pos[i].y, pos[i].z);
    fprintf(out, "END\n");
    fclose(out);
  }

  // Compute the current in nanoamperes.
  double current(double t) const {
    double curr = 0.0;
    double dt = timeLast - t;

    for (int i = 0; i < num; i++) {
      Vector3 d = wrap(pos[i]-posLast[i]);
      curr += part[type[i]].charge*d.z/(sysDim.z*dt)*1.60217733e-1;
    }
    return curr;
  }

  // Compute the current in nanoamperes for a restricted segment (-segZ < z < segZ).
  double currentSegment(double t, double segZ) const {
    double curr = 0.0;
    double dt = timeLast - t;

    for (int i = 0; i < num; i++) {
      double z0 = posLast[i].z;
      double z1 = pos[i].z;

      // Ignore carriers outside the range for both times.
      if (fabs(z0) > segZ && fabs(z1) > segZ) continue;

      // Cut the pieces outside the range.
      if (z0 < -segZ) z0 = -segZ;
      if (z1 < -segZ) z1 = -segZ;
      if (z0 > segZ) z0 = segZ;
      if (z1 > segZ) z1 = segZ;
      
      double dz = wrap(z1 - z0, sysDim.z);
      curr += part[type[i]].charge*dz/(sysDim.z*dt)*1.60217733e-1;
    }
    return curr;
  }

  void newCurrent() const {
    FILE* out = fopen(outCurrFile, "w");
    fclose(out);
  }

  void writeCurrent(double t) const {
    FILE* out = fopen(outCurrFile, "a");
    fprintf(out, "%.10g %.10g %d\n", 0.5*(t-timeLast), current(t), num);
    fclose(out);
  }

  void writeCurrentSegment(double t, double segZ) const {
    FILE* out = fopen(outCurrFile, "a");
    fprintf(out, "%.10g %.10g %d\n", 0.5*(t-timeLast), currentSegment(t, segZ), num);
    fclose(out);
  }

  static String makePdbLine(const String& tempLine, int index, const String& segName, int resId, 
		     const String& name, Vector3 r, double beta) {
    char s[128];

    String record("ATOM  ");
    sprintf(s, "     %5i ", index);
    String si = String(s).range(-6,-1);
    if (name.length() == 4) sprintf(s, "%s   ", name.val());
    else sprintf(s, " %s   ", name.val());
    String nam = String(s).range(0,3);
    String temp0 = tempLine.range(16,21);
  
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

    sprintf(s, "%s    ", segName.val());
    String seg = String(s).range(0,3);

    String ret(record);
    ret.add(si);
    ret.add(nam);
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

  void newPdbTraj(const char* fileName) {
    char s[128];

    sprintf(s, "CRYST1   %.3f   %.3f   %.3f  90.00  90.00  90.00 P 1           1\n", sysDim.x, sysDim.y, sysDim.z);
    String sysLine(s);

    sprintf(s, "REMARK   name %s\nREMARK   frameTime %.10g ns\n", outputName.val(), outputPeriod*timestep);
    String remarkLine(s);
    
    String tempLine("ATOM      1  CA  MOL S   1      -6.210  -9.711   3.288  0.00  0.00      ION");
    String line;

    FILE* out = fopen(fileName, "w");
    fprintf(out, "%s", sysLine.val());
    fprintf(out, "%s", remarkLine.val());

    for (int i = 0; i < num; i++) {
      line = makePdbLine(tempLine, i, part[type[i]].name, i, part[type[i]].name, pos[i], 0.0);
      fprintf(out, "%s",  line.val());
      fprintf(out, "\n");
    }
    fprintf(out, "END\n");
    fclose(out);
  }

  void writePdbTraj(const char* fileName) {
    String tempLine("ATOM      1  CA  MOL S   1      -6.210  -9.711   3.288  0.00  0.00      ION");
    String line;

    FILE* out = fopen(fileName, "a");
    for (int i = 0; i < num; i++) {
      line = makePdbLine(tempLine, i, part[type[i]].name, i, part[type[i]].name, pos[i], 0.0);
      fprintf(out, "%s", line.val());
      fprintf(out, "\n");
    }
    fprintf(out, "END\n");
    fclose(out);
  }


private:
  // Output variables
  char outCurrFile[256];
  char outTrajFile[256];
  char outPdbFile[256];

  // Integrator variables
  Vector3 sysDim;
  Random randoGen;
  ComputeForce* internal;
  Vector3* forceInternal;

  // Particle variables
  int numCap; // max number of particles
  int num; // current number of particles
  int num0; // mean number of particles
  Vector3* pos; //  position of each particle
  Vector3* pos1; // next position of each particle 
  int* type; // type of each particle
  String* name; // name of each particle
  Vector3* posLast; // used for current computation
  double timeLast; // used with posLast
  
  // Global parameters
  String outputName;
  double timestep;
  long int steps;
  long int seed;
  String inputCoordinates;
  int numberFluct;
  int interparticleForce;
  int tabulatedPotential;
  int fullElect;
  double kT;
  double coulombConst;
  double electricField;
  double cutoff;
  int outputPeriod;
  int outputPdb;
  double currentSegmentZ;
  int numberFluctPeriod;

  // Other parameters.
  int decompPeriod;
  double switchLen;
  double switchStart;
  double maxInitialPot;
  int numCapFactor;
  double initialZ;

  // Particle parameters.
  BrownianParticleType* part;
  int numParts;
  String* partGridFile;
  String* partTableFile;
  int* partTableIndex0;
  int* partTableIndex1;

  void remember(double t) {
    timeLast = t;
    for (int i = 0; i < num; i++) posLast[i] = pos[i];
  }
  
  void deleteParticles(IndexList& p) {
    int n = 0;
    for (int i = 0; i < num; i++) {
      pos[n] = pos[i];
      pos1[n] = pos1[i];
      type[n] = type[i];
      if (p.find(i) == -1) n++;
    }
    num = n;
  }

  void addParticles(int n, int typ) {
    if (num + n > numCap) n = numCap - num;

    for (int i = num; i < num + n; i++) {
      pos[i] = findPos(typ, initialZ);
      pos1[i] = pos[i];
      type[i] = typ;
    }
    num += n;
  }

  void numberFluctuation() {
    for (int p = 0; p < numParts; p++) {
      // Count the current number of this particle.
      int numberCurr = 0;
      IndexList ind;
      for (int i = 0; i < num; i++)
	if (type[i] == p) {
	  numberCurr++;
	  ind.add(i);
	}
      
      // Determine the new number for this particle.
      int n = randoGen.poisson(part[p].num);
      //printf("%d %d %d\n", p, numberCurr, n);

      // If the number is the same nothing needs to be done.
      if (n == numberCurr) continue;

      if (n < numberCurr) {
	// We need to delete particles.
	// Delete the farthest.
	IndexList farth(farthest(ind, numberCurr - n));
	deleteParticles(farth);
      } else {
	// We need to add particles.
	addParticles(n - numberCurr, p);
      }
    }
  }

  // Read the parameter file.
  int readParameters(const char* configFile) {
    Reader config(configFile);
    printf("Read config file %s.\n", configFile);
    
    const int numParams = config.length();
    numParts = config.countParameter("particle");
    part = new BrownianParticleType[numParts];
    partGridFile = new String[numParts];
    partTableFile = new String[numParts*numParts];
    partTableIndex0 = new int[numParts*numParts];
    partTableIndex1 = new int[numParts*numParts];

    int currPart = -1;
    int currTab = -1;
    for (int i = 0; i < numParams; i++) {
      String param = config.getParameter(i);
      String value = config.getValue(i);

      if (param == String("outputName") || param==String("outName")) outputName = value;
      else if (param == String("timestep")) timestep = strtod(value.val(), NULL);
      else if (param == String("steps")) steps = atol(value.val());
      else if (param == String("seed")) seed = atoi(value.val());
      else if (param == String("inputCoordinates")) inputCoordinates = value.val();
      else if (param == String("numberFluct")) numberFluct = atoi(value.val());
      else if (param == String("numberFluctPeriod")) numberFluctPeriod = atoi(value.val());
      else if (param == String("interparticleForce")) interparticleForce = atoi(value.val());
      else if (param == String("fullElect")) fullElect = atoi(value.val());
      else if (param == String("tabulatedPotential")) tabulatedPotential = atoi(value.val());
      else if (param == String("kT")) kT = strtod(value.val(), NULL);
      else if (param == String("coulombConst")) coulombConst = strtod(value.val(), NULL);
      else if (param == String("electricField")) electricField = strtod(value.val(), NULL);
      else if (param == String("cutoff")) cutoff = strtod(value.val(), NULL);
      else if (param == String("outputPeriod")) outputPeriod = atoi(value.val());
      else if (param == String("outputPdb")) outputPdb = atoi(value.val());
      else if (param == String("currentSegmentZ")) currentSegmentZ = strtod(value.val(), NULL);
      else if (param == String("particle")) { currPart++; part[currPart] = BrownianParticleType(value); }
      else if (param == String("num")) part[currPart].num = atoi(value.val());
      else if (param == String("gridFile")) partGridFile[currPart] = value;
      else if (param == String("diffusion")) part[currPart].diffusion = strtod(value.val(), NULL);
      else if (param == String("charge")) part[currPart].charge = strtod(value.val(), NULL);
      else if (param == String("radius")) part[currPart].radius = strtod(value.val(), NULL);
      else if (param == String("eps")) part[currPart].eps = strtod(value.val(), NULL);
      else if (param == String("tabulatedFile")) {currTab++; readTableFile(value, currTab);}
      else printf("WARNING: Unrecognized keyword `%s'.\n", param.val());
    }
   
    // Write the parameters.
    printf("\nParameters: \n");
    printf("outputName %s\n", outputName.val());
    printf("timestep %.10g\n", timestep);
    printf("steps %ld\n", steps);
    printf("seed %ld\n", seed);
    printf("inputCoordinates %s\n", inputCoordinates.val());
    printf("numberFluct %d\n", numberFluct);
    printf("numberFluctPeriod %d\n", numberFluctPeriod);
    printf("interparticleForce %d\n", interparticleForce);
    printf("fullElect %d\n", fullElect);
    printf("tabulatedPotential %d\n", tabulatedPotential);
    printf("kT %.10g\n", kT);
    printf("coulombConst %.10g\n", coulombConst);
    printf("electricField %.10g\n", electricField);
    printf("cutoff %.10g\n", cutoff);
    printf("outputPeriod %d\n", outputPeriod);
    printf("outputPdb %d\n", outputPdb);
    printf("currentSegmentZ %.10g\n", currentSegmentZ);

    // Write the particles.
    printf("\nParticles:\n");
  
    for (int i = 0; i < numParts; i++) {
      printf("particle %s\n", part[i].name.val());
      printf("num %d\n", part[i].num);
      printf("gridFile %s\n", partGridFile[i].val());
      printf("diffusion %.10g\n", part[i].diffusion);
      printf("charge %.10g\n", part[i].charge);
      printf("radius %.10g\n", part[i].radius);
      printf("eps %.10g\n\n", part[i].eps);
    }
    return numParams;
  }
};
