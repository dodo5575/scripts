///////////////////////////////////////////////////////////////////////
// Brownian dynamics base class
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef BROWNTOWN_H
#define BROWNTOWN_H

#include "Grid.H"

class  BrownTown {
public:
  BrownTown(const Grid& pot0, const Grid& gate0, double diffusion0, double kT0, double dt0) : pot(pot0), gate(gate0), diffusion(diffusion0), kT(kT0), dt(dt0), reflectX(1), reflectY(1), reflectZ(1), magicForce(0.0) {
  }

  void setReflection(int refX, int refY, int refZ) {
    reflectX = refX;
    reflectY = refY;
    reflectZ = refZ;
  }
  
  void setPeriodic(int perX, int perY, int perZ) {
    reflectX = !perX;
    reflectY = !perY;
    reflectZ = !perZ;
  }
  
  Vector3 step(Vector3 r0, Vector3 extForce, Vector3 rando, bool gated) {
    Vector3 f;

    if (gated) f = gate.interpolateForce(r0);
    else f = pot.interpolateForce(r0);
    
    Vector3 r = r0 + (f+extForce)*diffusion*dt/kT + sqrt(2*diffusion*dt)*rando;
    //printf("%s\n", rando.toString().val());
    magicForce = f;

    // Apply boundary conditions.
    Vector3 o = pot.getOrigin();
    Vector3 d = pot.getDestination();
    Vector3 l = pot.getExtent();

    if (reflectX) {
      // Reflecting boundaries.
      if (r.x < o.x) r.x = 2*o.x - r.x;
      if (r.x >= d.x) r.x = 2*d.x - r.x;
    } else {
      // Periodic boundaries.
      if (r.x < o.x) r.x += l.x;
      if (r.x >= d.x) r.x -= l.x;
    }
    
    if (reflectY) {
      // Reflecting boundaries.
      if (r.y < o.y) r.y = 2*o.y - r.y;
      if (r.y >= d.y) r.y = 2*d.y - r.y;
    } else {
      // Periodic boundaries.
      if (r.y < o.y) r.y += l.y;
      if (r.y >= d.y) r.y -= l.y;
    }
    
    if (reflectZ) {
      // Reflecting boundaries.
      if (r.z < o.z) r.z = 2*o.z - r.z;
      if (r.z >= d.z) r.z = 2*d.z - r.z;
    } else {
      // Periodic boundaries.
      if (r.z < o.z) r.z += l.z;
      if (r.z >= d.z) r.z -= l.z;
    }
    return r;
  }

public:
  Vector3 magicForce;

 private:
  double diffusion;
  Grid pot;
  Grid gate;
  double kT;
  double dt;
  int reflectX, reflectY, reflectZ;
};
#endif
