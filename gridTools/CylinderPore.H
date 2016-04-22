// Author: Jeff Comer <jcomer2@illinois.edu>

#include "useful.H"

#ifndef CYLINDERPORE_H
#define CYLINDERPORE_H

class CylinderPore {
public:
  CylinderPore(double poreLength, double poreDiamX, double poreDiamY, double cornerRad0) {
    endZ = 0.5*poreLength;
    cornerRad = cornerRad0;
    cornerZ = 0.5*poreLength - cornerRad;
    cornerRad2 = cornerRad*cornerRad;

    sx = 0.5 * poreDiamX;
    sy = 0.5 * poreDiamY;
    sx2 = sx*sx;
    sy2 = sy*sy;
  }

  bool in(Vector3 r) const {
    if (fabs(r.z) > endZ) return false;

    if (fabs(r.z) < cornerZ) {
      return (r.x*r.x/sx2 + r.y*r.y/sy2 > 1.0);
    }

    double dz = fabs(r.z) - cornerZ;
    double cornSx = sx + cornerRad - sqrt(cornerRad2 - dz*dz);
    double cornSy = sy + cornerRad - sqrt(cornerRad2 - dz*dz);
 
    return (r.x*r.x/(cornSx*cornSx) + r.y*r.y/(cornSy*cornSy) > 1.0);
  }

  bool onSurface(Vector3 r) {
    if (fabs(r.z) > endZ - 2*cornerRad) return true;
    return false;
  }

  double surfaceWeight(Vector3 r) {
    const double z = fabs(r.z);
    const double z0 = endZ - 2*cornerRad;
    if (z < z0) return 0.0;
    if (z > endZ) return 1.0;
    return (z0 - z)/(2.0*cornerRad);
  }

  double surfaceWeightSide(Vector3 r) {
    double cornSx = sx + cornerRad;
    double cornSy = sy + cornerRad;
    double ell = r.x*r.x/(cornSx*cornSx) + r.y*r.y/(cornSy*cornSy);
    
    if (ell <= 1.0) return 0.0;

    double w = 6.25*cornerRad*cornerRad/(sx*sy);
    if (ell >= 1.0 + w || cornerRad == 0.0) return 1.0;
    return (ell - 1.0)/w;
  }

private:
  double endZ, cornerZ, cornerRad, cornerRad2;
  double sx, sy, sx2, sy2;
};
#endif
