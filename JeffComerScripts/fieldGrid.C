#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "useful.H"
#include "Scatter.H"
#include "Grid.H"

using namespace std;

class FieldGrid {
public:
  FieldGrid(double eFieldMax0, double radXY0, double radZ0, Vector3 origin0)
    : eFieldMax(eFieldMax0), radXY(radXY0), radZ(radZ0), origin(origin0) {    
  }

  double gauss(Vector3 r) {
    Vector3 d = r - origin;
    double s2 = (d.x*d.x + d.y*d.y)/(radXY*radXY);
    double z2 = d.z*d.z/(radZ*radZ);
    return eFieldMax*exp(-0.5*(s2+z2));
  }

  double polyDrop(Vector3 r) {
    Vector3 d = r - origin;
    if (d.z <= -radZ) return 0.0;
    if (d.z >= radZ) return -4.0/3.0*eFieldMax*radZ;
    
    double x = d.z/radZ;
    return eFieldMax*radZ/3.0*(-2.0 - 3.0*x + x*x*x);
  }
  
  double biasDrop(Vector3 r) {
    Vector3 d = r - origin;
    if (d.z <= -radZ) return 0.0;
    if (d.z >= radZ) return -eFieldMax;
    
    double x = d.z/radZ;
    return 0.25*eFieldMax*(-2.0 - 3.0*x + x*x*x);
  }

  double drop(Vector3 r) {
    const double pi = 4.0*atan(1.0);
    const double b = 9.4;
    const double phi0 = pi*2.0*radZ*eFieldMax/b;
    
    Vector3 d = r - origin;

    return phi0/pi*atan(-b*d.z/(2.0*radZ));
  }

  double triangle(Vector3 r) {
    Vector3 d = r - origin;
    double s = sqrt(d.x*d.x + d.y*d.y)/radXY;
    double z = d.z/radZ;
    
    if (fabs(s) > radXY) return 0.0;
    if (fabs(d.z) > radZ) return 0.0;
    
  }

  double potential(Vector3 r) {
    Vector3 d = r - origin;
    if (fabs(d.z) > radZ) return 0.0;
    double s2 = (d.x*d.x + d.y*d.y)/(radXY*radXY);
    if (s2 > 1.0) return 0.0;
    double z = d.z/radZ;
    double z2 = z*z;
    
    double phi = eFieldMax;
    phi *= 1.0 - 2.0*s2 + s2*s2;
    //phi *= radZ*(8.0/15.0 + z2 - 2.0/3.0*z*z2 + 1.0/5.0*z2*z2*z);
    
    return phi;
  }
 
private:
  double eFieldMax;
  double radXY;
  double radZ;
  Vector3 origin;
};

///////////////////////////////////////////////////////////////////////
// Driver
int mainGrid(int argc, char* argv[]) {
  if ( argc != 9 ) {
    printf("Usage: %s basisFile bias radXY radZ ox oy oz outFile\n", argv[0]);
    printf("You entered %d parameters.\n", argc-1);
    return 0;
  }
  const double res = 0.1;

  // Load the basis vectors.
  Scatter cell(argv[1]);
  Matrix3 basis(cell.topMatrix());
  
  // Get the parameters.
  double eFieldMax = strtod(argv[2],NULL);
  double radXY = strtod(argv[3],NULL);
  double radZ = strtod(argv[4],NULL);
  double ox = strtod(argv[5],NULL);
  double oy = strtod(argv[6],NULL);
  double oz = strtod(argv[7],NULL);
  Vector3 origin = Vector3(ox,oy,oz);

  // Make the geometry object.
  FieldGrid geo(eFieldMax, radXY, radZ, origin);
  double dl = (radZ < radXY)?(res*radZ):(res*radXY);
  
  // Create Grid.
  Grid g(basis, dl);
 
  // Make the mask.
  const int n = g.length();
  for (int i = 0; i < n; i++) {
    Vector3 r = g.getPosition(i);
    g.setValue(i, geo.biasDrop(r));
    //g.setValue(i, geo.gauss(r));
  }

  char comments[256];
  sprintf(comments, "fieldGrid:  eFieldMax %g, radZ %g, radXY %g", eFieldMax, radZ, radXY);
  printf("%s\n", comments);
  g.write(argv[argc-1], comments);
  
  return 0;
}

int main(int argc, char* argv[]) {
  return mainGrid(argc, argv);
}
