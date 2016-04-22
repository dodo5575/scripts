///////////////////////////////////////////////////////////////////////
// Brownian dynamics base class
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef COMPUTEFORCE_H
#define COMPUTEFORCE_H

#include "Grid.H"
#include "BrownianParticleType.H"
#include "CellDecomposition.H"
#include "TabulatedPotential.H"

class ComputeForce {
public:
  ComputeForce(int num0, const BrownianParticleType* part, int numParts0, const Grid& g, double switchStart0, double switchLen0, double electricConst0) : 
    num(num0),
    numParts(numParts0),
    sys(g.getExtent()), switchStart(switchStart0),
    switchLen(switchLen0), electricConst(electricConst0),
    decomp(g.getBox(), g.getOrigin(), switchStart0+switchLen0) {

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
    tablePot = new TabulatedPotential*[numParts*numParts];
    for (int i = 0; i < numParts*numParts; i++) tablePot[i] = NULL;
    
    // Make the cell decomposition.
    neigh = new IndexList[num];
  }

  ~ComputeForce() {
    for (int i = 0; i < numParts; i++) {
      delete[] tableEps[i];
      delete[] tableRad6[i];
      delete[] tableAlpha[i];
    }
    delete[] tableEps;
    delete[] tableRad6;
    delete[] tableAlpha;

    for (int j = 0; j < numParts*numParts; j++) {
      if (tablePot[j] != NULL) {
	delete tablePot[j];
	tablePot[j] = NULL;
      }
    }
    delete[] tablePot;

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

  bool addTabulatedPotential(String fileName, int type0, int type1) {
    if (type0 < 0 || type0 >= numParts) return false;
    if (type1 < 0 || type1 >= numParts) return false;

    int ind = type0 + type1*numParts;
    int ind1 = type1 + type0*numParts;

    if (tablePot[ind] != NULL) {
      delete tablePot[ind];
      tablePot[ind] = NULL;
    }
    if (tablePot[ind1] != NULL) delete tablePot[ind1];
    
    tablePot[ind] = new TabulatedPotential(fileName);
    tablePot[ind1] = new TabulatedPotential(*tablePot[ind]);

    return true;
  }

  void decompose(const Vector3* pos) {
    // Reset the cell decomposition.
    decomp.clearCells();
    decomp.decompose(pos, num);

    // Regenerate the neighbor lists.
    for (int i = 0; i < num; i++) neigh[i] = decomp.neighborhood(pos[i]);
  }

  IndexList decompDim() const {
    IndexList ret;
    ret.add(decomp.getNx());
    ret.add(decomp.getNy());
    ret.add(decomp.getNz());
    return ret;
  }

  double decompCutoff() const {
    return decomp.getCutoff();
  }

  double computeFull(Vector3* force, const Vector3* pos, const int* type) const {
    double energy = 0.0;
    // Zero the force.
    for (int i = 0; i < num; i++) force[i] = Vector3(0.0);
    
    // Compute the force for all pairs.
    for (int i = 0; i < num-1; i++) {
      for (int j = i + 1; j < num; j++) {
	double alpha = tableAlpha[type[i]][type[j]];
	double eps = tableEps[type[i]][type[j]];
	double rad6 = tableRad6[type[i]][type[j]];
	Vector3 dr = wrap(pos[j] - pos[i]);

	EnergyForce fc = coulombForceFull(dr, alpha);
	EnergyForce fh = softcoreForce(dr, eps, rad6);

	force[i] += fc.f + fh.f;
	force[j] -= fc.f + fh.f;
	energy += fc.e + fh.e;
	//printf("TAB: %.10g %.10g\n", dr.length(), dr.dot(fc.f)/dr.length());
	//printf("TAB: %.10g %.10g\n", dr.length(), fc.e+fh.e);
      }
    }

    return energy;
  }

  double computeSoftcoreFull(Vector3* force, const Vector3* pos, const int* type) const {
    double energy = 0.0;
    // Zero the force.
    for (int i = 0; i < num; i++) force[i] = Vector3(0.0);
    
    // Compute the force for all pairs.
    for (int i = 0; i < num-1; i++) {
      for (int j = i + 1; j < num; j++) {
	double eps = tableEps[type[i]][type[j]];
	double rad6 = tableRad6[type[i]][type[j]];
	Vector3 dr = wrap(pos[j] - pos[i]);

	EnergyForce fh = softcoreForce(dr, eps, rad6);
	
	force[i] += fh.f;
	force[j] -= fh.f;
	energy += fh.e;
	//printf("TAB: %.10g %.10g\n", dr.length(), dr.dot(fh.f)/dr.length());
      }
    }

    return energy;
  }

  double compute(Vector3* force, const Vector3* pos, const int* type) const {
    double energy = 0.0;
    
    for (int i = 0; i < num; i++) {
      // Zero the force.
      force[i] = Vector3(0.0);

      // Loop through the neighbors.
      for (int n = 0; n < neigh[i].length(); n++) {
	int j = neigh[i].get(n);
	if (j == i) continue;
	double alpha = tableAlpha[type[i]][type[j]];
	double eps = tableEps[type[i]][type[j]];
	double rad6 = tableRad6[type[i]][type[j]];
	Vector3 dr = wrap(pos[j] - pos[i]);

	EnergyForce fc = coulombForce(dr, alpha, switchStart, switchLen);
	EnergyForce fh = softcoreForce(dr, eps, rad6);
	
	force[i] += fc.f + fh.f;
	energy += 0.5*(fc.e + fh.e);
	//printf("TAB: %.10g %.10g\n", dr.length(), dr.dot(fc.f + fh.f)/dr.length());
	//printf("TAB: %.10g %.10g\n", dr.length(), fc.e + fh.e);
      }
    }

    return energy;
  }

  double computeSoftcore(Vector3* force, const Vector3* pos, const int* type) const {
    double energy = 0.0;
    for (int i = 0; i < num; i++) {
      // Zero the force.
      force[i] = Vector3(0.0);
      
      // Loop through the neighbors.
      for (int n = 0; n < neigh[i].length(); n++) {
	int j = neigh[i].get(n);
	if (j == i) continue;
	double eps = tableEps[type[i]][type[j]];
	double rad6 = tableRad6[type[i]][type[j]];
	Vector3 dr = wrap(pos[j] - pos[i]);

	EnergyForce fh = softcoreForce(dr, eps, rad6);   
	
	force[i] += fh.f;
	energy += 0.5*fh.e;
      }
    }
    return energy;
  }

  double computeTabulated(Vector3* force, const Vector3* pos, const int* type) {
    double energy = 0.0;
    for (int i = 0; i < num; i++) {
      // Zero the force.
      force[i] = Vector3(0.0);
      
      // Loop through the neighbors.
      for (int n = 0; n < neigh[i].length(); n++) {
	int j = neigh[i].get(n);
	if (j == i) continue;
	int ind = type[i] + type[j]*numParts;
	if (tablePot[ind] == NULL) continue;
	Vector3 dr = wrap(pos[j] - pos[i]);
	
	EnergyForce ft = tablePot[ind]->compute(dr);

	force[i] += ft.f;
	energy += ft.e;
      }
    }
    return energy;
  }

  double computeTabulatedFull(Vector3* force, const Vector3* pos, const int* type) {
    energy = 0.0;
    // Zero the force.
    for (int i = 0; i < num; i++) force[i] = Vector3(0.0);
    
    // Compute the force for all pairs.
    for (int i = 0; i < num-1; i++) {
      for (int j = i + 1; j < num; j++) {
	int ind = type[i] + type[j]*numParts;
	if (tablePot[ind] == NULL) continue;
	Vector3 dr = wrap(pos[j] - pos[i]);

	EnergyForce ft = tablePot[ind]->compute(dr);
	//printf("TAB: %.10g %.10g\n", dr.length(), ft.e);

	force[i] += ft.f;
	force[j] -= ft.f;
	energy += ft.e;
      }
    }
    return energy;
  }

  static EnergyForce coulombForce(Vector3 r, double alpha, double start, double len) {
    double d = r.length();

    if (d >= start + len) return EnergyForce();
    if (d <= start) {
      double energy = alpha/d - alpha/start + 0.5*alpha/(start*start)*len;
      Vector3 force = -alpha/(d*d*d)*r;
      return EnergyForce(energy, force);
    }

    // Switching.
    double c = alpha/(start*start);
    double energy = 0.5/len*c*(start + len - d)*(start + len - d);
    Vector3 force = -c*(1.0 - (d - start)/len)/d*r;

    return EnergyForce(energy, force);
  }

  static EnergyForce coulombForceFull(Vector3 r, double alpha) {
    double d = r.length();
    
    return EnergyForce(alpha/d, -alpha/(d*d*d)*r);
  }

  static EnergyForce softcoreForce(Vector3 r, double eps, double rad6) {
    const double d2 = r.length2();
    const double d6 = d2*d2*d2;
  
    if (d6 < rad6) 
      return EnergyForce(eps*(rad6*rad6/(d6*d6) - 2.0*rad6/d6) + eps, 
			 -12*eps*(rad6*rad6/(d6*d6*d2) - rad6/(d6*d2))*r);
    return EnergyForce();
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
  int num;
  int numParts;
  double** tableEps;
  double** tableRad6;
  double** tableAlpha;
  Vector3 sys;
  IndexList* neigh;
  double switchStart, switchLen, electricConst;
  CellDecomposition decomp;
  TabulatedPotential** tablePot;
  int numTablePots;
  double energy;
};
#endif
