///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef TABULATEDPOTENTIAL_H
#define TABULATEDPOTENTIAL_H

#include "useful.H"
#define SWITCH_COUNT 20

class EnergyForce {
public:
  EnergyForce(double energy, Vector3 force) : e(energy), f(force) {
  }

  EnergyForce(): e(0.0), f(0.0) {
  }

  double e;
  Vector3 f;
};

class TabulatedPotential {
public:
  TabulatedPotential() {
    n = 2;
    r0 = 0.0;
    dr = 1.0;
    r1 = r0 + n*dr;
    
    v0 = new double[n];
    v0[0] = 0.0;
    v0[1] = 0.0;
    interpolate();
  }

  TabulatedPotential(const char* fileName) {
    FILE* inp = fopen(fileName, "r");
    char line[256];

    const int numLines = countValueLines(fileName);
    double* r = new double[numLines];
    double* v = new double[numLines];

    int count = 0;
    while (fgets(line, 256, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
      
      String s(line);
      int numTokens = s.tokenCount();
      if (numTokens != 2) {
	printf("Warning: Invalid tabulated potential file line: %s\n", line);
	continue;
      }

      String* tokenList = new String[numTokens];
      s.tokenize(tokenList);
      if (tokenList == NULL) {
	printf("Warning: Invalid tabulated potential file line: %s\n", line);
	continue;
      }
      r[count] = strtod(tokenList[0], NULL);
      v[count] = strtod(tokenList[1], NULL);
      count++;

      delete[] tokenList;
    }

    init(r, v, count);
    interpolate();
    delete[] r;
    delete[] v;
  }

  TabulatedPotential(const TabulatedPotential& tab) {
    n = tab.n;
    dr = tab.dr;
    r0 = tab.r0;
    r1 = tab.r1;

    v0 = new double[n];
    v1 = new double[n];
    v2 = new double[n];
    v3 = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = tab.v0[i];
      v1[i] = tab.v1[i];
      v2[i] = tab.v2[i];
      v3[i] = tab.v3[i];
    }
  }

  TabulatedPotential(const double* dist, const double* pot, int n0) {
    init(dist, pot, n0);
    interpolate();
  }

  ~TabulatedPotential() {
    delete[] v0;
    delete[] v1;
    delete[] v2;
    delete[] v3;
  }

  static int countValueLines(const char* fileName) {
    FILE* inp = fopen(fileName, "r");
    char line[256];
    int count = 0;

    while (fgets(line, 256, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
      
      count++;
    }
    fclose(inp);

    return count;
  }

  EnergyForce compute(Vector3 r) {
    double d = r.length();
    Vector3 rUnit = -r/d;
    int home = int(floor((d - r0)/dr));

    // Everything is zero after dr*(n+SWITCH_COUNT).
    if (home >= n + SWITCH_COUNT) return EnergyForce();

    // Switch the energy and force to zero.
    if (home >= n) {
      const int last = n-1;
      const double s = SWITCH_COUNT*dr;
      const double d0 = n*dr;
      const double e0 = v3[last] + v2[last] + v1[last] + v0[last];
      const double f0 = -(3.0*v3[last] + 2.0*v2[last] + v1[last])/dr;
          
      double x = (d - d0)/s;
      double x2 = x*x;
      double f0s = f0*s;
      double energy = e0 - f0s*x + (3.0*e0-2*f0s)*x2 + (-2.0*e0+f0s)*(x*x2);
      Vector3 force = (f0 - 2.0*(3.0*e0-2*f0s)*x/s - 3.0*(-2.0*e0+f0s)*x2/s)*rUnit;
      return EnergyForce(energy,force);
    }
        
    // Maintain a constant force below the first point in the table.
    if (home < 0) {
      const double e0 = v0[0];
      const double f0 = -v1[0]/dr;

      double energy = -f0*(d - r0) + e0;
      Vector3 force = f0*rUnit;
      return EnergyForce(energy,force);
    }

    double homeR = home*dr + r0;
    double w = (d - homeR)/dr;
   
    // Interpolate.
    double energy = v3[home]*w*w*w + v2[home]*w*w + v1[home]*w + v0[home];
    Vector3 force = -(3.0*v3[home]*w*w + 2.0*v2[home]*w + v1[home])*rUnit/dr;
    printf("FORCE: %.10g %.10g\n", d, force.dot(rUnit));
    return EnergyForce(energy,force);
  }

private:
  double* pot;
  double* v0;
  double* v1;
  double* v2;
  double* v3;
  int n;
  double dr;
  double r0, r1;
 
  void init(const double* dist, const double* pot, int n0) {
    n = abs(n0);
    dr = dist[1]-dist[0];
    r0 = dist[0];
    r1 = r0 + n*dr;

    v0 = new double[n];
    for (int i = 0; i < n; i++) v0[i] = pot[i];
  }
  
  void interpolate() {
    v1 = new double[n];
    v2 = new double[n];
    v3 = new double[n];

    for (int i = 0; i < n; i++) {
      int i0 = i - 1;
      int i1 = i;
      int i2 = i + 1;
      int i3 = i + 2;

      if (i0 < 0) i0 = 0;
      if (i2 >= n) i2 = n-1;
      if (i3 >= n) i3 = n-1;

      v3[i] = 0.5*(-v0[i0] + 3*v0[i1] - 3*v0[i2] + v0[i3]);
      v2[i] = 0.5*(2*v0[i0] - 5*v0[i1] + 4*v0[i2] - v0[i3]);
      v1[i] = 0.5*(-v0[i0] + v0[i2]);
    }
  }
 
};
#endif
