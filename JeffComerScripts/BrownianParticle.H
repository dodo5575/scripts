// Configuration file reader
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef BROWNIANPARTICLE_H
#define BROWNIANPARTICLE_H

class BrownianParticle {
public:
  BrownianParticle() : type(0), r(0.0) {
  }

public:
  int type;
  double r;
};

#endif 
