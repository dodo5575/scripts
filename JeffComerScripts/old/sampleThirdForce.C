/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995-2007 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

class Vector3;
class Matrix3;
class IndexList;
class ThirdForce;
class Grid;

// class Vector3
// Operations on 3D double vectors
//
class Vector3 {
public:
  Vector3() {}

  Vector3(double x0, double y0, double z0) {
    x = x0;
    y = y0;
    z = z0;
  }

  const Vector3 operator-() const {
    Vector3 v;
    v.x = -x;
    v.y = -y;
    v.z = -z;
    return v;
  }

  const Vector3 operator+(const Vector3 w ) const {
    Vector3 v;
    v.x = x + w.x;
    v.y = y + w.y;
    v.z = z + w.z;
    return v;
  }

  const Vector3 operator-(const Vector3 w ) const {
    Vector3 v;
    v.x = x - w.x;
    v.y = y - w.y;
    v.z = z - w.z;
    return v;
  }
  
  const Vector3 operator*(double s) const {
    Vector3 v;
    v.x = s*x;
    v.y = s*y;
    v.z = s*z;
    return v;
  }

  double length() const { return sqrt(x*x + y*y + z*z); }
  double length2() const { return x*x + y*y + z*z; }
  void toString(char* s) const {
    sprintf(s, "%.10g %.10g %.10g", x, y, z);
  }
  double x,y,z;
};

Vector3 operator*(double s, Vector3 v) {
  v.x *= s;
  v.y *= s;
  v.z *= s;
  return v;
}

Vector3 operator/(Vector3 v, double s) {
  v.x /= s;
  v.y /= s;
  v.z /= s;
  return v;
}

// class Matrix3
// Operations on 3D double matrices
//
class Matrix3 {
public:
  Matrix3() {}

  Matrix3(Vector3 ex, Vector3 ey, Vector3 ez) {
    exx = ex.x;
    eyx = ex.y;
    ezx = ex.z;
    exy = ey.x;
    eyy = ey.y;
    ezy = ey.z;
    exz = ez.x;
    eyz = ez.y;
    ezz = ez.z;
    
  }

  const Matrix3 operator*(double s) const {
    Matrix3 m;
    m.exx = s*exx;
    m.exy = s*exy;
    m.exz = s*exz;
    m.eyx = s*eyx;
    m.eyy = s*eyy;
    m.eyz = s*eyz;
    m.ezx = s*ezx;
    m.ezy = s*ezy;
    m.ezz = s*ezz;

    return m;
  }

  const Matrix3 operator-() const {
    Matrix3 m;
    m.exx = -exx;
    m.exy = -exy;
    m.exz = -exz;
    m.eyx = -eyx;
    m.eyy = -eyy;
    m.eyz = -eyz;
    m.ezx = -ezx;
    m.ezy = -ezy;
    m.ezz = -ezz;

    return m;
  }

  Matrix3 transpose() const {
    Matrix3 m;
    m.exx = exx;
    m.exy = eyx;
    m.exz = ezx;
    m.eyx = exy;
    m.eyy = eyy;
    m.eyz = ezy;
    m.ezx = exz;
    m.ezy = eyz;
    m.ezz = ezz;

    return m;
  }

  Matrix3 inverse() const {
    Matrix3 m;
    double det = (exx*(eyy*ezz-eyz*ezy) - exy*(eyx*ezz-eyz*ezx) + exz*(eyx*ezy-eyy*ezx));

    m.exx = (eyy*ezz - eyz*ezy)/det;
    m.exy = -(exy*ezz - exz*ezy)/det;
    m.exz = (exy*eyz - exz*eyy)/det;
    m.eyx = -(eyx*ezz - eyz*ezx)/det;
    m.eyy = (exx*ezz - exz*ezx)/det;
    m.eyz = -(exx*eyz - exz*eyx)/det;
    m.ezx = (eyx*ezy - eyy*ezx)/det;
    m.ezy = -(exx*ezy - exy*ezx)/det;
    m.ezz = (exx*eyy - exy*eyx)/det;

    return m;
  }

 Vector3 transform(Vector3 v) const {
   Vector3 w;
   w.x = exx*v.x + exy*v.y + exz*v.z;
   w.y = eyx*v.x + eyy*v.y + eyz*v.z;
   w.z = ezx*v.x + ezy*v.y + ezz*v.z;
   return w;
  }
  
  Vector3 a() const {return Vector3(exx,eyx,ezx);}
  Vector3 b() const {return Vector3(exy,eyy,ezy);}
  Vector3 c() const {return Vector3(exz,eyz,ezz);}

  void toString(char* s) const {
    sprintf(s, "%2.8f %2.8f %2.8f\n%2.8f %2.8f %2.8f\n%2.8f %2.8f %2.8f",
		   exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz);
  }

  double exx, exy, exz;
  double eyx, eyy, eyz;
  double ezx, ezy, ezz;
};

Matrix3 operator*(double s, Matrix3 m) { 
  m.exx *= s;
  m.exy *= s;
  m.exz *= s;
  m.eyx *= s;
  m.eyy *= s;
  m.eyz *= s;
  m.ezx *= s;
  m.ezy *= s;
  m.ezz *= s;
  return m;
}

Matrix3 operator/(Matrix3 m, double s) {
  m.exx /= s;
  m.exy /= s;
  m.exz /= s;
  m.eyx /= s;
  m.eyy /= s;
  m.eyz /= s;
  m.ezx /= s;
  m.ezy /= s;
  m.ezz /= s;
  return m;
}

// class IndexList
// A growable list of integers, for holding indices of atoms
//
class IndexList {
public:
  IndexList() {
    num=0;
    maxnum=16;
    list = new int[16];
  }
  
  ~IndexList() {
    if (list != 0)
      delete [] list;
  }
  
  void add(const int val) {
    // If we need more space, allocate a new block that is 1.5 times larger
    // and copy everything over
    if (num == maxnum) {
      int oldmax = maxnum;
      maxnum = (maxnum*3)/2 + 1;
      int* oldlist = list;
      list = new int[maxnum];
      int i;
      for(i=0;i<num;i++) {
        list[i] = oldlist[i];
      }
      delete [] oldlist;
    }
    
    // We should have enough space now, add the value
    list[num] = val;
    num++;
  }

  int size() {
    return num;
  }
  
  int get(const int i) const { return list[i]; }

private:
  int num;
  int maxnum;
  int* list;
};

class Grid {
public:
  void write(const char* fileName, const char* comments) {
    // Open the file.
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s\n.",fileName);
      exit(-1);
    }

    int size = nx*ny*nz;

    // Write the header.
    fprintf(out, "# %s\n", comments);
    fprintf(out, "object 1 class gridpositions counts %d %d %d\n", nx, ny, nz);
    fprintf(out, "origin %.10g %.10g %.10g\n", origin.x, origin.y, origin.z);
    fprintf(out, "delta %.10g %.10g %.10g\n", delta.exx, delta.eyx, delta.ezx);
    fprintf(out, "delta %.10g %.10g %.10g\n", delta.exy, delta.eyy, delta.ezy);
    fprintf(out, "delta %.10g %.10g %.10g\n", delta.exz, delta.eyz, delta.ezz);
    fprintf(out, "object 2 class gridconnections counts %d %d %d\n", nx, ny, nz);
    fprintf(out, "object 3 class array type double rank 0 items %d data follows\n", size);
    
    // Write the data.
    int penultima = 3*(size/3);
    int mod = size - penultima;

    int i;
    for (i = 0; i < penultima; i+=3) {
      fprintf(out, "%.10f %.10f %.10f\n", pot[i], pot[i+1], pot[i+2]);
    }
    if (mod == 1) {
      fprintf(out, "%.10f\n", pot[size-1]);
    } else if (mod == 2) {
      fprintf(out, "%.10f %.10f\n", pot[size-2], pot[size-1]);
    }
    fclose(out);
  }


  Vector3 origin;
  Matrix3 delta;
  int nx, ny, nz;
  double* pot;
};

class ThirdForceGrid {
public:
  ThirdForceGrid(const char* coordFile, const char* basisFile) {
    loaded = false;
    // Read the basis vectors.
    printf("Loading the basis from %s...\n", basisFile);
    Vector3 b[3];
    readCoordinates(basisFile, 3, b);
    basis.exx = b[0].x;
    basis.eyx = b[0].y;
    basis.ezx = b[0].z;
    basis.exy = b[1].x;
    basis.eyy = b[1].y;
    basis.ezy = b[1].z;
    basis.exz = b[2].x;
    basis.eyz = b[2].y;
    basis.ezz = b[2].z;
    char s[128];
    basis.toString(s);
    printf("\nBasis:\n%s\n\n", s);

    // Define the periodic basis vectors.
    period[0] = -basis.a();
    period[1] = basis.a();
    period[2] = -basis.b();
    period[3] = basis.b();
    period[4] = -basis.c();
    period[5] = basis.c();

    // Count the number of atoms.
    nAtoms = countCoordinates(coordFile);
    printf("Reading %d atoms from %s...\n", nAtoms, coordFile);

    // Read the coordinates.
    loaded = true;
    pos = new Vector3[nAtoms];
    readCoordinates(coordFile, nAtoms, pos);
    printf("Done.\n");
  }

  ~ThirdForceGrid() {
    // Deallocate the array if necessary.
    if (loaded) {
      delete[] pos;
    }
  }

  int countCoordinates(const char* fileName) {
    int nRead = 3;
    int n = -1;
    double x, y, z;

    // Open the file.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      printf("Couldn't open file %s\n.",fileName);
      exit(-1);
    }

    // Count the lines.
    while (nRead == 3) {
      n++;
      nRead = fscanf(inp,"%lf %lf %lf",&x,&y,&z);
    }
    
    fclose(inp);
    return n;
  }

  double computePotential(double r) {
    if (r < radius0) return radius0 - r + 0.5*sigma;
    if (r > radius0 + sigma) return 0.0;
    
    double dr = r - radius0;
    return dr*(0.5*dr/sigma - 1.0) + 0.5*sigma;
  }

  double potential(Vector3 r) {
    int j, p;
    double pot = 0.0;
    double dSqMin, dSq;
    Vector3 d;

    for (j = 0; j < nAtoms; j++) {
      // Get the identity distance.
      d = r - pos[j];
      dSqMin = d.length2();

      // Find the minimum distance among the images.
      for (p = 0; p < 6; p++) {
	d = r + period[p] - pos[j];
	dSq = d.length2();
	if (dSq < dSqMin) dSqMin = dSq;
      }
     
      pot += computePotential(sqrt(dSqMin));
    }

    return pot;
  }

  Grid samplePotential(double rad0, double sig, double dx) {
    char s[128];
    radius0 = rad0;
    sigma = sig;

    int na = int(ceil(basis.a().length()/dx));
    int nb = int(ceil(basis.b().length()/dx));
    int nc = int(ceil(basis.c().length()/dx));
        
    Vector3 a = basis.a()/na;
    Vector3 b = basis.b()/nb;
    Vector3 c = basis.c()/nc;
    Vector3 diag = basis.a() + basis.b() + basis.c();

    Grid g;
    g.origin = diag*-0.5;
    g.nx = na;
    g.ny = nb;
    g.nz = nc;
    g.delta = Matrix3(a,b,c);
    int size = g.nx*g.ny*g.nz;
    
    printf("\nSampling potential with resolution of %.10g.\n", dx);
    printf("radius0 = %.10g, sigma = %.10g.\n", radius0, sigma);
    printf("Grid size is %d.\n", size);
    printf("Grid dimensions are %d %d %d.\n", g.nx, g.ny, g.nz);
    g.origin.toString(s);
    printf("Grid origin is %s.\n", s);
    g.delta.toString(s);
    printf("Grid basis:\n%s\n", s);

    // Allocate the potential array.
    g.pot = new double[size];

    // Sample.
    printf("\nSampling the potential...\n");
    int ia, ib, ic, j;
    Vector3 r;
    double pot;
    int n = 0;
    for (ia = 0; ia < na; ia++) {
      for (ib = 0; ib < nb; ib++) {
	for (ic = 0; ic < nc; ic++) {
	  r = g.origin + a*ia + b*ib + c*ic;
	  
	  g.pot[n] = potential(r);
	  n++;
	}
      }

      if (ia % 5 == 0) printf("%2.1f percent complete\n", (100.0*n)/size);
    }

    return g;
  }

double radius0, sigma;
private:
  int nAtoms;
  bool loaded;
  Vector3* pos; 
  Matrix3 basis;
  Vector3 period[6];

  void readCoordinates(const char* fileName, int n, Vector3* r) {
    // Open the file.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      printf("Couldn't open file %s\n.",fileName);
      exit(-1);
    }
    
    int i;
    int nRead;
    double x, y, z;
    for (i = 0; i < n; i++) {
      nRead = fscanf(inp,"%lf %lf %lf",&x,&y,&z);
      if (nRead != 3) {
	printf("Error reading atom %d (%d)\n.",i,nRead);
	exit(-2);
      }

      r[i].x = x;
      r[i].y = y;
      r[i].z = z; 
    }

    fclose(inp);
  }
};


int main(int argc,char *argv[])
{
  if ( argc != 6 ) {
    printf(
      "Usage: %s coordFile basisFile radius0 sigma outFile\n",
      argv[0]);
    return 0;
  }

  // Extract the parameters.
  double radius0 = strtod(argv[3],NULL);
  double sigma = strtod(argv[4],NULL);
  printf("\n********\nThird force grid: \n");
  printf("Creating a grid with radius0 = %.10g and sigma = %.10g...\n",
	 radius0, sigma);

  // Organize the file comments.
  char comments[128];
  sprintf(comments, "ThirdForce grid for %s with radius0 = %.10g and sigma = %.10g", argv[1], radius0, sigma);
  
  // Generate the grid.
  ThirdForceGrid third(argv[1], argv[2]);
  Grid g = third.samplePotential(radius0, sigma, 0.8);  
  printf("Writing the grid...\n");
  g.write(argv[5], comments);
  printf("Wrote %s successfully.\n", argv[5]);

  // Don't forget to deallocate.
  delete[] g.pot;
  return 0;
}



