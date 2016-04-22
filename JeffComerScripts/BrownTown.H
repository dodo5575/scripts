///////////////////////////////////////////////////////////////////////
// Brownian dynamics base class
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef BROWNTOWN_H
#define BROWNTOWN_H

#include "Grid.H"

class  BrownTown {
public:
  BrownTown(double kT0, double dt0, const Grid& g) : kT(kT0), dt(dt0), reflectX(1), reflectY(1), reflectZ(1), sys(g.getExtent()) {
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

  double wrap(double x, double l) {
    int image = int(floor(x/l));
    x -= image*l;
    
    if (x > 0.5*l) x -= l;
    return x;
  }

  Vector3 wrap(Vector3 r) {
    r.x = wrap(r.x, sys.x);
    r.y = wrap(r.y, sys.y);
    r.z = wrap(r.z, sys.z);
    return r;
  }

  Vector3 stepPeriodic(Vector3 r0, Vector3 force, Vector3 rando, double diffusion) {
    Vector3 r = r0 + force*diffusion*dt/kT + sqrt(2*diffusion*dt)*rando;
    // Apply boundary conditions.
    r = wrap(r);
    
    return r;
  }

  Vector3 step(Vector3 r0, Vector3 o, Vector3 d, Vector3 force, Vector3 rando, double diffusion) {
    Vector3 r = r0 + force*diffusion*dt/kT + sqrt(2*diffusion*dt)*rando;

    // Apply boundary conditions.
    const Vector3 l = d - o;

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

private:
  double kT;
  double dt;
  int reflectX, reflectY, reflectZ;
Vector3 sys;
};
#endif
