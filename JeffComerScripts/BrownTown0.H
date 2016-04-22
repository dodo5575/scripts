///////////////////////////////////////////////////////////////////////
// Brownian dynamics base class
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef BROWNTOWN_H
#define BROWNTOWN_H

#include "Grid1.H"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


class  BrownTown {
public:
  BrownTown(const Grid& pot0, double diffusion0, double kT0, double dt0, double z0) : pot(pot0), diffusion(diffusion0), kT(kT0), dt(dt0), edgeZ(z0), reflectX(1), reflectY(1), reflectZ(1), magicForce(0.0) {
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
  
  static Vector3 normalDeviate() {
    // Box-Mueller transformation
    double x1 = 2.0*(double(rand())/RAND_MAX - 0.5);
    double x2 = 2.0*(double(rand())/RAND_MAX - 0.5);
    
    double z1 = sqrt(-2*log(x1))*cos(2.0*M_PI*x2);
    double z2 = sqrt(-2*log(x1))*sin(2.0*M_PI*x2);
    
    double x3 = 2.0*(double(rand())/RAND_MAX - 0.5);
    double x4 = 2.0*(double(rand())/RAND_MAX - 0.5);
    
    double z3 = sqrt(-2*log(x1))*cos(2.0*M_PI*x2);

    return Vector3(z1,z2,z3);    
  }

  Vector3 step(Vector3 r0) {
  
  // select random number generator 
  r = gsl_rng_alloc (gsl_rng_mt19937);

	double sigma = 10; 
  for (int i = 0; i < 10; i++) 
    {
      double u = gsl_rng_uniform (r);
      double v = gsl_ran_gaussian(r, sigma); 
    }

  gsl_rng_free (r);



    //Vector3 rando = normalDeviate();
    Vector3 rando = gsl_ran_gaussian_ziggurat(const gsl_rng* r, 1.0);
    
    Vector3 f;
    if (r0.z < edgeZ) {
      f = pot.interpolateForce(r0);
    } else {
      f = Vector3(0.0);
    }
  
    Vector3 r = r0 + f*diffusion*dt/kT + sqrt(2*diffusion*dt)*rando;
    printf("%g\n", rando.toString().val());
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
  double kT;
  double dt;
  int reflectX, reflectY, reflectZ;
  double edgeZ;
};
#endif
