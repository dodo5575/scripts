// Configuration file reader
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef BROWNIANPARTICLETYPE_H
#define BROWNIANPARTICLETYPE_H

class BrownianParticleType {
public:
  BrownianParticleType() : name(), diffusion(1.0), charge(1.0), radius(1.0), eps(1.0), grid(NULL) {
  }
  BrownianParticleType(const String& name0) : name(name0), diffusion(1.0), charge(1.0), radius(1.0), eps(1.0), grid(NULL) {
  }

  ~BrownianParticleType() {
    if (grid != NULL) delete grid;
  }

  String name;
  int num;
  double diffusion;
  double charge;
  double radius;
  double eps;
  Grid* grid;
};

#endif 
