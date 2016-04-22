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

// A mod like wrapping function
double wrap(double n, double d) {
  return d - n*int(floor(d/n));
}

// class Vector3
// Operations on 3D double vectors
//
class Vector3 {
public:
  Vector3() {}
  
  Vector3(const Vector3& v) {
    x = v.x;
    y = v.y;
    z = v.z;
  }
  
  Vector3(double x0, double y0, double z0) {
    x = x0;
    y = y0;
    z = z0;
  }

  const Vector3& operator=(const Vector3& v) {
    x = v.x;
    y = v.y;
    z = v.z;
    return *this;
  }

  const Vector3 operator-() const {
    Vector3 v;
    v.x = -x;
    v.y = -y;
    v.z = -z;
    return v;
  }

  const Vector3 operator+(const Vector3& w ) const {
    Vector3 v;
    v.x = x + w.x;
    v.y = y + w.y;
    v.z = z + w.z;
    return v;
  }

  const Vector3 operator-(const Vector3& w ) const {
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

  Matrix3(double s) {
    exx = s;
    exy = 0;
    exz = 0;
    eyx = 0;
    eyy = s;
    eyz = 0;
    ezx = 0;
    ezy = 0;
    ezz = s;
  }

  Matrix3(const Vector3& ex, const Vector3& ey, const Vector3& ez) {
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

 Vector3 transform(const Vector3& v) const {
   Vector3 w;
   w.x = exx*v.x + exy*v.y + exz*v.z;
   w.y = eyx*v.x + eyy*v.y + eyz*v.z;
   w.z = ezx*v.x + ezy*v.y + ezz*v.z;
   return w;
  }
  
  Vector3 ex() const {return Vector3(exx,eyx,ezx);}
  Vector3 ey() const {return Vector3(exy,eyy,ezy);}
  Vector3 ez() const {return Vector3(exz,eyz,ezz);}

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
    list = new int[maxnum];
  }

  IndexList(const IndexList& l) {
    num = l.num;
    maxnum = num + 16;
    list = new int[maxnum];

    for(int i = 0; i < num; i++) list[i] = l.list[i];
  }
  
  ~IndexList() {
    delete[] list;
  }
  
  void add(const int val) {
    // If we need more space, allocate a new block that is 1.5 times larger
    // and copy everything over
    if (num == maxnum) {
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

  void operator=(const IndexList& l) {
    delete[] list;

    num = l.num;
    maxnum = num + 16;
    list = new int[maxnum];

    for(int i = 0; i < num; i++) list[i] = l.list[i];
  }

  void add(const IndexList& l) {
    int oldnum = num;
    num = num + l.num;

    if (num > maxnum) {
      maxnum = (num*3)/2 + 1;
      int* oldlist = list;
      list = new int[maxnum];
      
      for(int i = 0; i < oldnum; i++) {
	list[i] = oldlist[i];
      }
      delete[] oldlist;
    }

    for(int i = oldnum; i < num; i++) {
      list[i] = l.list[i];
    }
  }

  int length() {
    return num;
  }
  
  int get(const int i) const {return list[i];}
  
  void clear() {
    num=0;
    maxnum=16;
    delete[] list;
    list = new int[maxnum];
  }

private:
  int num;
  int maxnum;
  int* list;
};


///////////////////////////////////////////////////////////////////////
// Interesting stuff.

///////////////////////////////////////////////////////////////////////
// An array of positions.
class Scatter {
public:
  Scatter(const char* coordFile) {
    // Count the number of points.
    n = countCoordinates(coordFile);
    
    // Load the coordinates.
    r = new Vector3[n];
    readCoordinates(coordFile, n, r);
  }
  
  ~Scatter() {
    delete[] r;
  }

  Matrix3 topMatrix() const {
    if (n < 3) return Matrix3(1.0);
    return Matrix3(r[0], r[1], r[2]);
  }
  Vector3 get(int i) const {
    return r[i];
  }
  int length() const {
    return n;
  }

  Vector3 minBound() const {
    Vector3 ret = r[0];
    for (int i = 1; i < n; i++) {
      if (r[i].x < ret.x) ret.x = r[i].x;
      if (r[i].y < ret.y) ret.y = r[i].y;
      if (r[i].z < ret.z) ret.z = r[i].z;
    }
    return ret;
  }

  Vector3 maxBound() const {
    Vector3 ret = r[0];
    for (int i = 1; i < n; i++) {
      if (r[i].x > ret.x) ret.x = r[i].x;
      if (r[i].y > ret.y) ret.y = r[i].y;
      if (r[i].z > ret.z) ret.z = r[i].z;
    }
    return ret;
  }

  static int countCoordinates(const char* fileName) {
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

private:
  int n;
  Vector3* r;

  Scatter(const Scatter&){};

  // Read coordinates into a Vector array.
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


///////////////////////////////////////////////////////////////////////
// Grid base class.
class Grid {
public:
  Grid(Matrix3 basis0, Vector3 origin0, int nx0, int ny0, int nz0) {
    basis = basis0;
    origin = origin0;
    nx = abs(nx0);
    ny = abs(ny0);
    nz = abs(nz0);
    
    basisInv = basis.inverse();
    size = nx*ny*nz;
    val = new double[size];
    zero();
  }

  Grid(Matrix3 box, Vector3 origin0, double dx) {
    dx = fabs(dx);
    
    // Tile the grid into the system box.
    // The grid spacing is always a bit larger than dx.
    nx = int(floor(box.ex().length()/dx))-1;
    ny = int(floor(box.ey().length()/dx))-1;
    nz = int(floor(box.ez().length()/dx))-1;
    basis = Matrix3(box.ex()/nx, box.ey()/ny, box.ez()/nz);
    origin = origin0;

    basisInv = basis.inverse();
    size = nx*ny*nz;
    val = new double[size];
    zero();
  }

  virtual ~Grid() {
    delete[] val;
  }

  virtual void zero() {
    for (int i = 0; i < size; i++) {
      val[i] = 0.0;
    }
  }

  virtual void write(const char* fileName, const char* comments) {
    // Open the file.
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s\n.",fileName);
      exit(-1);
    }

    // Write the header.
    fprintf(out, "# %s\n", comments);
    fprintf(out, "object 1 class gridpositions counts %d %d %d\n", nx, ny, nz);
    fprintf(out, "origin %.10g %.10g %.10g\n", origin.x, origin.y, origin.z);
    fprintf(out, "delta %.10g %.10g %.10g\n", basis.exx, basis.eyx, basis.ezx);
    fprintf(out, "delta %.10g %.10g %.10g\n", basis.exy, basis.eyy, basis.ezy);
    fprintf(out, "delta %.10g %.10g %.10g\n", basis.exz, basis.eyz, basis.ezz);
    fprintf(out, "object 2 class gridconnections counts %d %d %d\n", nx, ny, nz);
    fprintf(out, "object 3 class array type double rank 0 items %d data follows\n", size);
    
    // Write the data.
    int penultima = 3*(size/3);
    int mod = size - penultima;

    int i;
    for (i = 0; i < penultima; i+=3) {
      fprintf(out, "%.10f %.10f %.10f\n", val[i], val[i+1], val[i+2]);
    }
    if (mod == 1) {
      fprintf(out, "%.10f\n", val[size-1]);
    } else if (mod == 2) {
      fprintf(out, "%.10f %.10f\n", val[size-2], val[size-1]);
    }
    fclose(out);
  }


protected:
  Vector3 origin;
  Matrix3 basis;
  int nx, ny, nz;
  int size;
  Matrix3 basisInv;
  double* val;

private:
  Grid(const Grid&) {}
};


///////////////////////////////////////////////////////////////////////
// Cell decomposition of points.
struct SpaceCell {
  IndexList point;
  IndexList neigh;
  int ix, iy, iz;
  int index;
};

class CellDecomposition : public Grid {
public:
  // Make an empty cell decomposition with the desired geometry.
  CellDecomposition(Matrix3 box0, Vector3 origin0, double cutoff) : Grid(box0, origin0, cutoff) {

    cut = cutoff;
    box = box0;
    
    cell = new SpaceCell[size];
    clearCells();
    populateNeighbors();
  }

  virtual ~CellDecomposition() {
    delete[] cell;
  }

  // Find the cell index for a given position.
  int getCell(const Vector3& r) {
    // Transform.
    Vector3 l = basisInv.transform(r-origin);
    // Wrap into the home box.
    int cx = int(floor(wrap(nx, l.x)));
    int cy = int(floor(wrap(ny, l.y)));
    int cz = int(floor(wrap(nz, l.z)));
    
    return cz + cy*nz + cx*nz*ny;
  }

  // Restart with an empty decomposition.
  void clearCells() {
    int j = 0;

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  cell[j].ix = ix;
	  cell[j].iy = iy;
	  cell[j].iz = iz;
	  cell[j].index = j;
	  cell[j].point.clear();
	  val[j] = 0.0;

	  j++;
	}
      }
    }
  }

  // Generate the neighbor lists.
  void populateNeighbors() {
    int j = 0;

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  cell[j].neigh = getNeighbors(j);
	  j++;
	}
      }
    }
  }

  // Fill the neighbor list with the cell's neighbors (including itself).
  IndexList getNeighbors(int j) {
    int ix, iy, iz;
    IndexList ret;

    for (ix = -1; ix <= 1; ix++) {
      for (iy = -1; iy <= 1; iy++) {
	for (iz = -1; iz <= 1; iz++) {
	  int cx = cell[j].ix + ix;
	  if (cx >= nx) cx -= nx;
	  if (cx < 0) cx += nx;

	  int cy = cell[j].iy + iy;
	  if (cy >= ny) cy -= ny;
	  if (cy < 0) cy += ny;

	  int cz = cell[j].iz + iz;
	  if (cz >= nz) cz -= nz;
	  if (cz < 0) cz += nz;

	  ret.add(cz + cy*nz + cx*nz*ny);
	}
      }
    }

    return ret;
  }

  void addPoint(const Vector3& r, int pointIndex) {
    int j = getCell(r);
    cell[j].point.add(pointIndex);
    val[j] += 1.0;
  }

  void decompose(const Scatter& pos) {
    const int n = pos.length();
    clearCells();

    for (int i = 0; i < n; i++) {
      addPoint(pos.get(i), i);
    }
  }

  int countPoints() {
    int num = 0;
    for (int j = 0; j < size; j++) {
      num += cell[j].point.length();
    }
    return num;
  }

protected:
  double cut; // cutoff length
  Matrix3 box; // periodic box of the system
  SpaceCell* cell; // contains cell attributes
};

///////////////////////////////////////////////////////////////////////
// 


///////////////////////////////////////////////////////////////////////
// Drivers
int cellDecomp(int argc,char *argv[])
{
  if ( argc != 5 ) {
    printf("Usage: %s coordFile basisFile cutoff outFile\n", argv[0]);
    return 0;
  }
  double cutoff = strtod(argv[3],NULL);
  char s[128];
  
  Scatter pos(argv[1]);
  printf("Pos length: %i\n", pos.length());
  
  Scatter sysVec(argv[2]);
  Matrix3 sys = sysVec.topMatrix();
  sys.toString(s);
  printf("System length: %i\n", sysVec.length());
  printf("System:\n%s\n\n", s);
  Vector3 origin = (sys.ex() + sys.ey() + sys.ez())*-0.5;
  
  // Generate the cell decomposition.
  CellDecomposition cell(sys, origin, cutoff);
  //Grid cell(sys, origin, cutoff);
  cell.decompose(pos);
  printf("CellDecomposition grid for %s with cutoff = %.10g\n", argv[1], cutoff);
  printf("CellDecomposition contains %i points\n", cell.countPoints());
  
  
  // Organize the file comments.
  char comments[128];
  sprintf(comments, "CellDecomposition grid for %s with cutoff = %.10g", argv[1], cutoff);
  cell.write(argv[4], comments);

  return 0;
}


int main(int argc,char *argv[])
{
  if ( argc != 7 ) {
    printf(
      "Usage: %s coordFile basisFile dx radius0 sigma outFile\n",
      argv[0]);
    return 0;
  }

  // Extract the parameters.
  double dx = strtod(argv[3],NULL);
  double radius0 = strtod(argv[4],NULL);
  double sigma = strtod(argv[5],NULL);
  printf("\n********\nThird force grid: \n");
  printf("The resolution is %.10g.\n", dx);
  printf("Creating a grid with radius0 = %.10g and sigma = %.10g...\n",
	 radius0, sigma);

  // Get the system dimensions.
  Scatter sysVec(argv[2]);
  Matrix3 sys = sysVec.topMatrix();
  sys.toString(s);
  printf("System dimensions:\n%s\n", s);
  Vector3 origin = (sys.ex() + sys.ey() + sys.ez())*-0.5;
  origin.toString(s);
  printf("System origin:%s\n", s);

  // Generate the cell decomposition.
  double cutoff = 2.0*(radius0+sigma);
  printf("\nCellDecomposition grid for %s with cutoff = %.10g\n", argv[1], cutoff);
  printf("Generating the cell decomposition...");
  CellDecomposition cell(sys, origin, cutoff);
  //Grid cell(sys, origin, cutoff);
  cell.decompose(pos);
  printf("CellDecomposition contains %i points\n", cell.countPoints());
  

  // Organize the file comments.
  char comments[128];
  sprintf(comments, "ThirdForce grid for %s with dx = %.10g, radius0 = %.10g, and sigma = %.10g", argv[1], dx, radius0, sigma);
  
  // Generate the grid.
  ThirdForceGrid third(argv[1], argv[2]);
  Grid g = third.samplePotential(radius0, sigma, dx);  
  printf("Writing the grid...\n");
  g.write(argv[6], comments);
  printf("Wrote %s successfully.\n", argv[6]);

  // Don't forget to deallocate.
  delete[] g.pot;
  return 0;
}




