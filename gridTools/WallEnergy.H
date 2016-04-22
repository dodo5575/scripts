// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef WALLENERGY_H
#define WALLENERGY_H

class WallEnergy {
public:
  // xStart is where the function begins to deviate from zero.
  // xMin, yMin are the position and value of the function minimum.
  // xRise, yRise are the position and value near the wall.
  WallEnergy(double xStart0, double xMin0, double yMin0, double xRise0, double yRise0) {
    xStart = xStart0;
    
    // Shift by the starting x value.
    double zMin = xMin0 - xStart;
    double yMin = yMin0;
    zRise = xRise0 - xStart;
    yRise = yRise0;

    // Compute the parameters.
    double zMin2 = zMin*zMin;
    double zRise2 = zRise*zRise;
    double denom = (zRise - zMin)*(zRise - zMin)*zRise2*zMin2;

    c = (zMin2*zMin2*yRise + 3*zRise2*zRise2*yMin - 4*zRise2*zRise*zMin*yMin)/denom;
    d = (4*yMin - 2*c*zMin2)/(zMin2*zMin);
    e = (c*zMin2 - 3*yMin)/(zMin2*zMin2);
    lastSlope = 2*c*zRise + 3*d*zRise2 + 4*e*zRise2*zRise;
  }

  // Actually compute the value.
  double energy(double x) {
    double z = x - xStart;
    
    if (z < 0) return 0.0;
    if (z > zRise) return lastSlope*(z - zRise) + yRise;
    
    double z2 = z*z;
    return c*z2 + d*z2*z + e*z2*z2;
  }

  double hardEnergy(double x) {
  double z = x - xStart;
    
    if (z < 0) return 0.0;
    return z*z*yRise/(zRise*zRise);
  }
  
private:
  double xStart, zRise, yRise, lastSlope;
  double c, d, e;
};
#endif
