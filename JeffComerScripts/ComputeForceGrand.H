///////////////////////////////////////////////////////////////////////
// Brownian dynamics base class
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef COMPUTEFORCEGRAND_H
#define COMPUTEFORCEGRAND_H

#include "Grid.H"
#include "BrownianParticle.H"
#include "CellDecomposition.H"

class  ComputeForceGrand {
public:
  ComputeForceGrand(int numCap0, const BrownianParticleType* part, int numParts0, const Grid& g, double switchStart0, double switchLen0, double electricConst0) : 
    numCap(numCap0),
    numParts(numParts0),
    sys(g.getExtent()), switchStart(switchStart0),
    switchLen(switchLen0), electricConst(electricConst0),
    decomp(g.getBox(), g.getOrigin(), 1.2*(switchStart0+switchLen0)) {

    // Allocate the parameter tables.
    tableEps = new double*[numParts];
    tableRad6 = new double*[numParts];
    tableAlpha = new double*[numParts];
    for (int i = 0; i < numParts; i++) {
      tableEps[i] = new double[numParts];
      tableRad6[i] = new double[numParts];
      tableAlpha[i] = new double[numParts];
    }    

    // Form the parameter tables.
    makeTables(part);

    // Make the cell decomposition.
    neigh = new IndexList[numCap];
  }

  ~ComputeForceGrand() {
    for (int i = 0; i < numParts; i++) {
      delete[] tableEps[i];
      delete[] tableRad6[i];
      delete[] tableAlpha[i];
    }
    delete[] tableEps;
    delete[] tableRad6;
    delete[] tableAlpha;

    delete[] neigh;
  }

  void makeTables(const BrownianParticleType* part) {
    for (int i = 0; i < numParts; i++) {
      for (int j = 0; j < numParts; j++) {
	tableEps[i][j] = sqrt(part[i].eps*part[j].eps);
	double r = part[i].radius + part[j].radius;
	tableRad6[i][j] = r*r*r*r*r*r;
	tableAlpha[i][j] = electricConst*part[i].charge*part[j].charge;
      }
    }
  }

  void decompose(const Vector3* pos, int num) {
    decomp.clearCells();
    decomp.decompose(pos, num);

    // Regenerate the neighbor lists.
    for (int i = 0; i < num; i++) neigh[i] = decomp.neighborhood(pos[i]);
  }

  void computeFull(Vector3* force, const Vector3* pos, const int* type, int num) const {
    // Zero the force.
    for (int i = 0; i < num; i++) force[i] = 0.0;
    
    // Compute the force for all pairs.
    for (int i = 0; i < num-1; i++) {
      for (int j = i + 1; j < num; j++) {
	double alpha = tableAlpha[type[i]][type[j]];
	double eps = tableEps[type[i]][type[j]];
	double rad6 = tableRad6[type[i]][type[j]];
	Vector3 dr = wrap(pos[j] - pos[i]);

	Vector3 fc = coulombForceFull(dr, alpha);
	Vector3 fh = hardcoreForce(dr, eps, rad6);

	force[i] += fc + fh;
	force[j] -= fc + fh;
      }
    }
  }

  void computeHardcoreFull(Vector3* force, const Vector3* pos, const int* type, int num) const {
    // Zero the force.
    for (int i = 0; i < num; i++) force[i] = 0.0;
    
    // Compute the force for all pairs.
    for (int i = 0; i < num-1; i++) {
      for (int j = i + 1; j < num; j++) {
	double alpha = tableAlpha[type[i]][type[j]];
	double eps = tableEps[type[i]][type[j]];
	double rad6 = tableRad6[type[i]][type[j]];
	Vector3 dr = wrap(pos[j] - pos[i]);

	Vector3 fh = hardcoreForce(dr, eps, rad6);
	force[i] += fh;
	force[j] -= fh;
      }
    }
  }

  void compute(Vector3* force, const Vector3* pos, const int* type, int num) const {
    for (int i = 0; i < num; i++) {
       // Zero the force.
      force[i] = 0.0;
      
      // Loop through the neighbors.
      for (int n = 0; n < neigh[i].length(); n++) {
	int j = neigh[i].get(n);
	if (j == i) continue;
	double alpha = tableAlpha[type[i]][type[j]];
	double eps = tableEps[type[i]][type[j]];
	double rad6 = tableRad6[type[i]][type[j]];
	Vector3 dr = wrap(pos[j] - pos[i]);

	Vector3 fc = coulombForce(dr, alpha, switchStart, switchLen);
	Vector3 fh = hardcoreForce(dr, eps, rad6);
      
	force[i] += fc + fh;
      }
    }
  }

  void computeHardcore(Vector3* force, const Vector3* pos, const int* type, int num) const {
    for (int i = 0; i < num; i++) {
       // Zero the force.
      force[i] = 0.0;
      
      // Loop through the neighbors.
      for (int n = 0; n < neigh[i].length(); n++) {
	int j = neigh[i].get(n);
	if (j == i) continue;
	double eps = tableEps[type[i]][type[j]];
	double rad6 = tableRad6[type[i]][type[j]];
	Vector3 dr = wrap(pos[j] - pos[i]);

	Vector3 fh = hardcoreForce(dr, eps, rad6);
      
	force[i] += fh;
      }
    }
  }

  static Vector3 coulombForce(Vector3 r, double alpha, double start, double len) {
    double d = r.length();

    if (d > start + len) return 0;
    if (d < start) return alpha/(d*d*d)*r;
    return -alpha/(d*d*d)*(1 - (d-start)/len)*r;
  }

  static Vector3 coulombForceFull(Vector3 r, double alpha) {
    double d = r.length();

    return -alpha/(d*d*d)*r;
  }

  static Vector3 hardcoreForce(Vector3 r, double eps, double rad6) {
    const double d = r.length();
    const double d2 = d*d;
    const double d6 = d2*d2*d2;
  
    if (d6 < rad6) return -12*eps*(rad6*rad6/(d6*d6*d2) - rad6/(d6*d2))*r;
    return 0.0;
  }
  
  static double wrap(double x, double l) {
    int image = int(floor(x/l));
    x -= image*l;
    
    if (x > 0.5*l) x -= l;
    return x;
  }

  Vector3 wrap(Vector3 r) const {
    r.x = wrap(r.x, sys.x);
    r.y = wrap(r.y, sys.y);
    r.z = wrap(r.z, sys.z);
    return r;
  }

private:
  int numCap;
  int numParts;
  double** tableEps;
  double** tableRad6;
  double** tableAlpha;
  Vector3 sys;
  IndexList* neigh;
  double switchStart, switchLen, electricConst;
  CellDecomposition decomp;
};
#endif
