// Use the weighted histogram analysis method.
// Author: Jeff Comer <jcomer2@illinois.edu>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define DEBUG

class Vector3;
class Matrix3;
class IndexList;

bool isReal(char c) {
  char num[] = "0123456789-eE.";
  
  for (int i = 0; i < 14; i++) {
    if (c == num[i]) return true;
  }
  return false;
}

int firstSpace(const char* s, int max) {
  for (int i = 0; i < max; i++) {
    if (s[i] == ' ') return i;
  }
  return -1;
}

// A mod like wrapping function
int wrap(int n, double d) {
  double q = (d - n*int(floor(d/n)));
  int w = int(floor(q));
  if (w < 0) w = n;
  if (w >= n) w = 0;
  return w;
}

// A classic growable string class.
class String {
public:
  String() {
    cap = 16;
    c = new char[cap];
    c[0] = '\0';
    len = 1;
  }

  String(const char* s) {
    len = strlen(s) + 1;
    cap = len;
    c = new char[cap];
    for (int i = 0; i < len; i++) c[i] = s[i];
  }

  String(const String& s) {
    len = s.len;
    cap = s.len;
    c = new char[cap];
    for (int i = 0; i < s.len; i++) c[i] = s.c[i];
  }

  String& operator=(const String& s) {
    delete[] c;

    len = s.len;
    cap = s.len;
    c = new char[cap];
    for (int i = 0; i < s.len; i++) c[i] = s.c[i];
    return *this;
  }

  String& operator=(const char* s) {
    delete[] c;
    len = strlen(s) + 1;
    cap = len;
    c = new char[cap];
    for (int i = 0; i < len; i++) c[i] = s[i];
    return *this;
  }

  void add(const char* s) {
    int n = strlen(s) + 1;
    len--;

    if (n + len > cap) grow(n + len);
    for (int i = 0; i < n; i++) c[i+len] = s[i];
    len += n;
  }

  void add(String& s) {
    len--;

    if (len + s.len > cap) grow(len + s.len);
    for (int i = 0; i < s.len; i++) c[i+len] = s[i];
    len += s.len;
  }
  
  int length() const {
    return len-1;
  }

  ~String() {
    delete[] c;
  }
  
  operator const char*() {
    return c;
  }

  const char* val() {
    return c;
  }
  
private:
  char* c;
  int cap, len;

  void grow(int n) {
    char* c0 = c;
    c = new char[n];
    cap = n;
    for (int i = 0; i < len; i++) c[i] = c0[i];
    delete[] c0;
  }
};

// class Vector3
// Operations on 3D double vectors
//
class Vector3 {
public:
  Vector3() {}

  Vector3(double s) {
    x = s;
    y = s;
    z = s;
  }

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

  static Vector3 random(double s) {
    Vector3 v;
    v.x = (double(rand())/RAND_MAX-0.5)*s;
    v.y = (double(rand())/RAND_MAX-0.5)*s;
    v.z = (double(rand())/RAND_MAX-0.5)*s;
    return v;
  }

  Vector3& operator=(const Vector3& v) {
    x = v.x;
    y = v.y;
    z = v.z;
    return *this;
  }

  Vector3& operator+=(const Vector3& v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }
  
  Vector3& operator-=(const Vector3& v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }
  
  Vector3& operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }
  
  Vector3& operator/=(double s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }

  Vector3 operator-() const {
    Vector3 v;
    v.x = -x;
    v.y = -y;
    v.z = -z;
    return v;
  }

  Vector3 operator+(const Vector3& w ) const {
    Vector3 v;
    v.x = x + w.x;
    v.y = y + w.y;
    v.z = z + w.z;
    return v;
  }

  Vector3 operator-(const Vector3& w ) const {
    Vector3 v;
    v.x = x - w.x;
    v.y = y - w.y;
    v.z = z - w.z;
    return v;
  }
  
  Vector3 operator*(double s) const {
    Vector3 v;
    v.x = s*x;
    v.y = s*y;
    v.z = s*z;
    return v;
  }

  double dot(const Vector3& w) const {
    return x*w.x + y*w.y + z*w.z;
  }

  Vector3 cross(const Vector3& w) const {
    Vector3 v;
    v.x = y*v.z - z*v.y;
    v.y = z*v.x - x*v.z;
    v.z = x*v.y - y*v.x;
    return v;
  }

  double length() const { return sqrt(x*x + y*y + z*z); }
  double length2() const { return x*x + y*y + z*z; }
  String toString() const {
    char s[128];
    sprintf(s, "%.10g %.10g %.10g", x, y, z);
    return String(s);
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

  Matrix3(double x, double y, double z) {
    exx = x;
    exy = 0;
    exz = 0;
    eyx = 0;
    eyy = y;
    eyz = 0;
    ezx = 0;
    ezy = 0;
    ezz = z;
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
    double det = exx*(eyy*ezz-eyz*ezy) - exy*(eyx*ezz-eyz*ezx) + exz*(eyx*ezy-eyy*ezx);

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

  double det() const {
    return exx*(eyy*ezz-eyz*ezy) - exy*(eyx*ezz-eyz*ezx) + exz*(eyx*ezy-eyy*ezx);
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

  String toString() const {
    char s[128];
    sprintf(s, "%2.8f %2.8f %2.8f\n%2.8f %2.8f %2.8f\n%2.8f %2.8f %2.8f",
		   exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz);
    return String(s);
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
    num = 0;
    maxnum = 16;
    lis = new int[maxnum];
  }

  IndexList(const IndexList& l) {
    num = l.num;
    maxnum = num + 16;
    lis = new int[maxnum];

    for(int i = 0; i < l.num; i++) lis[i] = l.lis[i];
  }
  
  ~IndexList() {
    delete[] lis;
  }
 
  void add(const int val) {
    // If we need more space, allocate a new block that is 1.5 times larger
    // and copy everything over
    if (num == maxnum) {
      maxnum = (maxnum*3)/2 + 1;
      int* oldlis = lis;
      lis = new int[maxnum];
      int i;
      for(i = 0; i < num; i++) {
        lis[i] = oldlis[i];
      }
      delete [] oldlis;
    }
    
    // We should have enough space now, add the value
    lis[num] = val;
    num++;
  }

  IndexList& operator=(const IndexList& l) {
    delete[] lis;

    num = l.num;
    maxnum = num + 16;
    lis = new int[maxnum];

    for(int i = 0; i < num; i++) lis[i] = l.lis[i];
    return *this;
  }

  void add(const IndexList& l) {
    int oldnum = num;
    num = num + l.num;

    if (num > maxnum) {
      maxnum = (num*3)/2 + 1;
      int* oldlis = lis;
      lis = new int[maxnum];
      
      for(int i = 0; i < oldnum; i++) lis[i] = oldlis[i];
      delete[] oldlis;
    }

    for(int i = 0; i < l.num; i++) lis[i+oldnum] = l.lis[i];
  }

  int length() const {
    return num;
  }
  
  int get(const int i) const {
#ifdef DEBUG 
    if (i < 0 || i >= num) {
      printf("Warning! IndexList::get out of bounds.\n");
      return 0;
    }
#endif
    return lis[i];
  }
  
  void clear() {
    num=0;
    maxnum=16;
    delete[] lis;
    lis = new int[maxnum];
  }

  String toString() const {
    String ret;
    char tmp[32];

    for (int i = 0; i < num; i++) {
      sprintf(tmp, "%i ", lis[i]);
      ret.add(tmp);
    }
    return ret;
  }

private:
  int num;
  int maxnum;
  int* lis;
};


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

  Scatter(const char* coordFile, double cutTime0, double cutTime1) {
    // Count the number of points.
    n = countTrajectory(coordFile, cutTime0, cutTime1);
    
    // Load the coordinates.
    r = new Vector3[n];
    readTrajectory(coordFile, n, r, cutTime0, cutTime1);
  }
  
  ~Scatter() {
    delete[] r;
  }

  Matrix3 topMatrix() const {
    if (n < 3) return Matrix3(1.0);
    return Matrix3(r[0], r[1], r[2]);
  }
  Vector3 get(int i) const {
#ifdef DEBUG 
    if (i < 0 || i >= n) {
      printf("Warning! Scatter::get out of bounds.\n");
      return Vector3(0.0);
    }
#endif
    return r[i];
  }
  virtual int length() const {
    return n;
  }

  virtual Vector3 minBound() const {
    Vector3 ret = r[0];
    for (int i = 1; i < n; i++) {
      if (r[i].x < ret.x) ret.x = r[i].x;
      if (r[i].y < ret.y) ret.y = r[i].y;
      if (r[i].z < ret.z) ret.z = r[i].z;
    }
    return ret;
  }

  virtual Vector3 maxBound() const {
    Vector3 ret = r[0];
    for (int i = 1; i < n; i++) {
      if (r[i].x > ret.x) ret.x = r[i].x;
      if (r[i].y > ret.y) ret.y = r[i].y;
      if (r[i].z > ret.z) ret.z = r[i].z;
    }
    return ret;
  }

  static int countCoordinates(const char* fileName) {
    int nRead;
    int n = 0;
    double x, y, z;
    char line[256];

   // Open the file.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
      exit(-1);
    }

    while (fgets(line, 256, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
    
      // Read values.
      nRead = sscanf(line, "%lf %lf %lf", &x, &y, &z);
      if (nRead >= 3) n++;
    }
    
    fclose(inp);
    return n;
  }

  static int countTrajectory(const char* fileName, double cutTime0, double cutTime1) {
    int nRead;
    int n = 0;
    double t, x, y, z;
    char line[256];

    // Open the file.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
      exit(-1);
    }

    while (fgets(line, 256, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
    
      // Read values.
      nRead = sscanf(line, "%lf %lf %lf %lf", &t, &x, &y, &z);
      if (nRead >= 4 && t >= cutTime0 && t< cutTime1) n++;
    }
    
    fclose(inp);
    return n;
  }

private:
  int n;
  Vector3* r;

  Scatter(const Scatter&){};

  // Read coordinates into a Vector array.
  virtual void readCoordinates(const char* fileName, int num, Vector3* r) {
    int nRead;
    int n = 0;
    double x, y, z;
    char line[256];

    // Open the file.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
      exit(-1);
    }

    while (fgets(line, 256, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
    
      // Read values.
      nRead = sscanf(line, "%lf %lf %lf %lf", &x, &y, &z);
      if (nRead >= 3) {
	r[n].x = x;
	r[n].y = y;
	r[n].z = z;
	n++;
	if (n >= num) break;
      }
    }
    
    fclose(inp);
  }

  // Read coordinates into a Vector array.
  virtual void readTrajectory(const char* fileName, int num, Vector3* r, double cutTime0, double cutTime1) {
    int nRead;
    int n = 0;
    double t, x, y, z;
    char line[256];

    // Open the file.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
      exit(-1);
    }

    while (fgets(line, 256, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
    
      // Read values.
      nRead = sscanf(line, "%lf %lf %lf %lf", &t, &x, &y, &z);
      if (nRead >= 4 && t >= cutTime0 && t < cutTime1) {
	r[n].x = x;
	r[n].y = y;
	r[n].z = z;
	n++;
	if (n >= num) break;
      }
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

  Grid(const Grid& g) {
    nx = g.nx;
    ny = g.ny;
    nz = g.nz;
    basis = g.basis;
    origin = g.origin;
    
    basisInv = g.basisInv;
    size = g.size;
    val = new double[size];
    for (int i = 0; i < size; i++) val[i] = g.val[i];
  }

  Grid& operator=(const Grid& g) {
    delete[] val;

    nx = g.nx;
    ny = g.ny;
    nz = g.nz;
    basis = g.basis;
    origin = g.origin;
    
    basisInv = g.basisInv;
    size = g.size;
    val = new double[size];
    for (int i = 0; i < size; i++) val[i] = g.val[i];

    return *this;
  }

  // Read a grid from a file.
  Grid(const char* fileName) {
    // Open the file.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      printf("Grid:Grid Couldn't open file %s\n.",fileName);
      exit(-1);
    }
    //printf("Reading dx file %s...\n", fileName);
    
    size = 0;
    nx = 0;
    ny = 0;
    nz = 0;
    basis = Matrix3(1.0);
    origin = Vector3(0.0);    

    int n = 0;
    double x, y, z;
    char line[256];
    int p, nRead;
    int deltaCount = 0;
    Vector3 base[3];
    while (fgets(line, 256, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
    
      if (isReal(line[0]) && n < size) {
	// Read grid values.
	nRead = sscanf(line, "%lf %lf %lf", &x, &y, &z);
	if (size > 0) {
	  switch(nRead) {
	  case 1:
	    val[n] = x;
	    n++;
	    break;
	  case 2:
	    val[n] = x;
	    val[n+1] = y;
	    n += 2;
	    break;
	  case 3:
	    val[n] = x;
	    val[n+1] = y;
	    val[n+2] = z;
	    n += 3;
	    break;
	  }
	}
      } else if (len > 5) {
	// Read the grid parameters.
	char start[6];
	for (int i = 0; i < 5; i++) {
	  start[i] = line[i];
	}

	if(strcmp("origi", start) == 0) {
	  // Get an origin line.
	  p = firstSpace(line, 256);
	  sscanf(&(line[p+1]), "%lf %lf %lf", &x, &y, &z);
	  origin = Vector3(x, y, z);
	  //printf("Origin: %.10g %.10g %.10g\n", x, y, z);
	} else if(strcmp("delta", start) == 0) {
	  // Get a delta matrix line.
	  p = firstSpace(line, 256);
	  sscanf(&(line[p+1]), "%lf %lf %lf", &x, &y, &z);
	  base[deltaCount] = Vector3(x, y, z);
	  //printf("Delta %d: %.10g %.10g %.10g\n", deltaCount, x, y, z);
	  if (deltaCount < 2) deltaCount = deltaCount + 1;
	} else if(strcmp("objec", start) == 0) {
	  // Get the system dimensions.
	  if (line[7] != '1') continue;
	  sscanf(line, "object 1 class gridpositions counts %d %d %d\n", &nx, &ny, &nz);
	  size = nx*ny*nz;
	  val = new double[size];
	  zero();
	  //printf("Size: %d %d %d\n", nx, ny, nz);
	}
      }
    }
    fclose(inp);

    basis = Matrix3(base[0], base[1], base[2]);
    basisInv = basis.inverse();
    if (size == 0) {
      printf("Grid:Grid Improperly formatted dx file %s.\n",fileName);
      exit(-1);
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
    fprintf(out, "origin %.14g %.14g %.14g\n", origin.x, origin.y, origin.z);
    fprintf(out, "delta %.14g %.14g %.14g\n", basis.exx, basis.eyx, basis.ezx);
    fprintf(out, "delta %.14g %.14g %.14g\n", basis.exy, basis.eyy, basis.ezy);
    fprintf(out, "delta %.14g %.14g %.14g\n", basis.exz, basis.eyz, basis.ezz);
    fprintf(out, "object 2 class gridconnections counts %d %d %d\n", nx, ny, nz);
    fprintf(out, "object 3 class array type double rank 0 items %d data follows\n", size);
    
    // Write the data.
    int penultima = 3*(size/3);
    int mod = size - penultima;

    int i;
    for (i = 0; i < penultima; i+=3) {
      fprintf(out, "%.14g %.14g %.14g\n", val[i], val[i+1], val[i+2]);
    }
    if (mod == 1) {
      fprintf(out, "%.14g\n", val[size-1]);
    } else if (mod == 2) {
      fprintf(out, "%.14g %.14g\n", val[size-2], val[size-1]);
    }
    fclose(out);
  }

  virtual ~Grid() {
    delete[] val;
  }

  virtual void zero() {
    for (int i = 0; i < size; i++) {
      val[i] = 0.0;
    }
  }
  
  virtual bool setValue(int j, double v) {
    if (j < 0 || j >= size) return false;
    val[j] = v;
    return true;
  }

  virtual bool setValue(int ix, int iy, int iz, double v) {
    if (ix < 0 || ix >= nx) return false;
    if (iy < 0 || iy >= ny) return false;
    if (iz < 0 || iz >= nz) return false;
    int j = iz + iy*nz + ix*ny*nz;

    val[j] = v;
    return true;
  }

  virtual double getValue(int j) const {
    if (j < 0 || j >= size) return 0.0;
    return val[j];
  }

  virtual double getValue(int ix, int iy, int iz) const {
    if (ix < 0 || ix >= nx) return 0.0;
    if (iy < 0 || iy >= ny) return 0.0;
    if (iz < 0 || iz >= nz) return 0.0;
    
    int j = iz + iy*nz + ix*ny*nz;
    return val[j];
  }

  virtual Vector3 getPosition(int j) const {
    int iz = j%nz;
    int iy = (j/nz)%ny;
    int ix = j/(nz*ny);

    return basis.transform(Vector3(ix, iy, iz)) + origin;
  }

  virtual int length() const {
    return size;
  }
  virtual Vector3 getOrigin() const {return origin;}
  virtual Matrix3 getBasis() const {return basis;}
  virtual int getNx() const {return nx;}
  virtual int getNy() const {return ny;}
  virtual int getNz() const {return nz;}
  virtual int getSize() const {return size;}
  virtual void setBasis(const Matrix3& b) {
    basis = b;
    basisInv = basis.inverse();
  }
  virtual void setOrigin(const Vector3& o) {
    origin = o;
  }
  virtual Vector3 getExtent() const {
    return basis.transform(Vector3(nx,ny,nz));
  }
  virtual double getCellVolume() const {
    return fabs(basis.det());
  }
  virtual Vector3 getVolume() const {
    return getCellVolume()*size;
  }

  virtual void shift(double s) {
    for (int i = 0; i < size; i++) val[i] += s;
  }

  virtual void scale(double s) {
    for (int i = 0; i < size; i++) val[i] *= s;
  }

  virtual void shiftToCenters() {
    origin += basis.transform(Vector3(0.5));
  }

  // Add a constant gradient to the grid.
  virtual void addGradient(Vector3 g) {
    int j = 0;

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  Vector3 dr = basis.transform(Vector3(ix,iy,iz));
	  double pot0 = dr.dot(g);
	  val[j] += pot0;
	  j++;
	}
      }
    }
  }

  // Get a profile along z, at a position defined by factorX and factorY.
  // factorX = factorY = 0 means that the profile starts from the origin.
  // factorX = factorY = 0.5 means the profile is along the center of the grid.
  virtual void profileZ(double factorX, double factorY, const char* fileName) {
    const int ix = int(floor((factorX*nx + 0.5)));
    const int iy = int(floor((factorY*ny + 0.5)));
    const int nynz = ny*nz;
    
    if (ix < 0 || ix >= nx || iy < 0 || iy >= ny) {
      printf("Error Grid::profileZ: Invalid x and y positions.");
      return;
    }
    
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s\n.",fileName);
      return;
    }

    for (int iz = 0; iz < nz; iz++) {
      int j = iz + iy*nz + ix*nynz;
      double v = val[j];
      double z = origin.z + iz*basis.ezz;
      fprintf(out, "%0.14g %0.14g\n", z, v);
    }
    fclose(out);
  }

  // Compute the average profile along z.
  virtual void averageProfileZ(const char* fileName) {
    const int nynz = ny*nz;
    
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s\n.",fileName);
      exit(-1);
    }

   
    for (int iz = 0; iz < nz; iz++) {
      double sum = 0.0;

      for (int ix = 0; ix < nx; ix++) {
	for (int iy = 0; iy < ny; iy++) {
	  int j = iz + iy*nz + ix*nynz;
	  sum += val[j];
	}
      }
      
      double v = sum/(nx*ny);
      double z = origin.z + iz*basis.ezz;
      fprintf(out, "%0.14g %0.14g\n", z, v);
    }

    fclose(out);
  }

  // Create profile along z averaged over a cylinder.
  virtual void cylinderProfileZ(const char* fileName, double radius) {
    const int nynz = ny*nz;
    const double rad2 = radius*radius;

    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s\n.",fileName);
      exit(-1);
    }

    for (int iz = 0; iz < nz; iz++) {
      double sum = 0.0;

      int n = 0;
      for (int ix = 0; ix < nx; ix++) {
	for (int iy = 0; iy < ny; iy++) {
	  Vector3 r = basis.transform(Vector3(ix,iy,iz)) + origin;

	  if (r.x*r.x + r.y*r.y < rad2) {
	    int j = iz + iy*nz + ix*nynz;
	    sum += val[j];
	    n++;
	  }
	}
      }
      
      double v = sum/n;
      double z = origin.z + iz*basis.ezz;
      fprintf(out, "%0.14g %0.14g\n", z, v);
    }

    fclose(out);
  }

  // Compute the average value over a box.
  virtual double averageRegion(Vector3 r0, Vector3 r1) {
    int i = 0;
    int n = 0;
    double sum = 0.0;

     for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  Vector3 r = basis.transform(Vector3(ix,iy,iz)) + origin;
	  if (r.x >= r0.x && r.x < r1.x && r.y >= r0.y && r.y < r1.y 
	      && r.z >= r0.z && r.z < r1.z) {
	    sum += val[i];
	    n++;
	  }
	  
	  i++;
	}
      }
     }

     return sum/n;
  }

  // Add two grids. The resulting grids has the dimensions of *this.
  // The added values come from the spatially nearest nodes from g.
  void addGrid(const Grid& g) {
    const int nynz = g.ny*g.nz;

    int i = 0;
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  Vector3 r = basis.transform(Vector3(ix,iy,iz)) + origin;
	  Vector3 rg = g.basisInv.transform(r-g.origin);
	  
	  // Find the nearest node in g.
	  int gx = int(floor(rg.x + 0.5));
	  int gy = int(floor(rg.y + 0.5));
	  int gz = int(floor(rg.z + 0.5));

	  // Check that this node lies in the domain of g.
	  bool good = true;
	  if (gx < 0) good = false;
	  if (gx >= g.nx) good = false;
	  if (gy < 0) good = false;
	  if (gy >= g.ny) good = false;
	  if (gz < 0) good = false;
	  if (gz >= g.nz) good = false;

	  // Add the value of the closest node.
	  if (good) {
	    int k = gz + gy*g.nz + gx*nynz;
	    val[i] += g.val[k];
	  }
	  i++;
	}
      }
    }
  }

  // Remove parts of the grid outside of the box defined by
  // r0 and r1.
  void crop(Vector3 r0, Vector3 r1) {
    const int nynz = ny*nz;

    // Find the start and end along x.
    int ix0 = -1;
    int ix1 = -1;
    for (int ix = 0; ix < nx; ix++) {
      Vector3 r = basis.transform(Vector3(ix,0,0)) + origin;
      if (ix0 < 0 && r.x >= r0.x) ix0 = ix;
      if (ix1 < 0 && r.x > r1.x) {
	ix1 = ix-1;
	break;
      }
    }
    if (ix0 < 0) ix0 = 0;
    if (ix1 < 0) ix1 = nx-1;

    // Find the start and end along y.
    int iy0 = -1;
    int iy1 = -1;
    for (int iy = 0; iy < ny; iy++) {
      Vector3 r = basis.transform(Vector3(0,iy,0)) + origin;
      if (iy0 < 0 && r.y >= r0.y) iy0 = iy;
      if (iy1 < 0 && r.y > r1.y) {
	iy1 = iy-1;
	break;
      }
    }
    if (iy0 < 0) iy0 = 0;
    if (iy1 < 0) iy1 = ny-1;

    // Find the start and end along z.
    int iz0 = -1;
    int iz1 = -1;
    for (int iz = 0; iz < nz; iz++) {
      Vector3 r = basis.transform(Vector3(0,0,iz)) + origin;
      if (iz0 < 0 && r.z >= r0.z) iz0 = iz;
      if (iz1 < 0 && r.z > r1.z) {
	iz1 = iz-1;
	break;
      }
    }
    if (iz0 < 0) iz0 = 0;
    if (iz1 < 0) iz1 = nz-1;

    int newNx = ix1-ix0+1;
    int newNy = iy1-iy0+1;
    int newNz = iz1-iz0+1;
    int newSize = newNx*newNy*newNz;
    double* v = new double[newSize];
    
    // Copy the appropriate data.
    int i = 0;
    for (int ix = ix0; ix <= ix1; ix++) {
      for (int iy = iy0; iy <= iy1; iy++) {
	for (int iz = iz0; iz <= iz1; iz++) {
	  int j = iz + iy*nz + ix*nynz;
	  v[i] = val[j];
	  i++;
	}
      }
    }

    // Determine the new origin and set the new members.
    origin += basis.transform(Vector3(ix0,iy0,iz0));
    nx = newNx;
    ny = newNy;
    nz = newNz;
    size = newSize;

    // Swap the pointers and deallocate the old array.
    double* valOld = val;
    val = v;
    delete[] valOld;
  }

  // Make a map of the minimum value along z.
  // Format is "x y z valMin".
  virtual void depthMapZ(const char* fileName) {
    FILE* out = fopen(fileName, "w");
    if (out == NULL) {
      printf("Grid:depthMapZ Couldn't open file %s\n.",fileName);
      exit(-1);
    }

    int i = 0;
    // Loop over x and y.
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	
	// Find the minimum value along z.
	int minIz = i;
	double minVal = val[i];
	for (int iz = 0; iz < nz; iz++) {
	  if (val[i] < minVal) {
	    minVal = val[i];
	    minIz = iz;
	  }
	  i++;
	}

	// Print the map of the minimum values.
	Vector3 r = basis.transform(Vector3(ix,iy,minIz)) + origin;
	fprintf(out, "%.14g %.14g %.14g %.14g\n", r.x, r.y, r.z, minVal);
      }
    }
    fclose(out);
  }

protected:
  Grid() {}
  Vector3 origin;
  Matrix3 basis;
  int nx, ny, nz;
  int size;
  Matrix3 basisInv;
  double* val;

private:
  //Grid(const Grid&) {}
  //Grid operator=(const Grid&) {}
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
  int getCell(const Vector3& r) const {
    // Transform.
    Vector3 l = basisInv.transform(r-origin);
    // Wrap into the home box.
    int cx = wrap(nx, l.x);
    int cy = wrap(ny, l.y);
    int cz = wrap(nz, l.z);
    
    int j = cz + cy*nz + cx*nz*ny;
#ifdef DEBUG 
    if (cx < 0 || cx >= nx) {
      printf("Warning! CellDecomposition::getCell cx out of bounds.\n");
      return 0;
    }
    if (cy < 0 || cy >= ny) {
      printf("Warning! CellDecomposition::getCell cy out of bounds.\n");
       return 0;
    }
    if (cz < 0 || cz >= nz) {
      printf("Warning! CellDecomposition::getCell cz out of bounds.\n");
      return 0;
    }
#endif
    return j;
  }
  
  IndexList getCellContents(int i) {
#ifdef DEBUG 
    if (i < 0 || i >= size) {
      printf("Warning! CellDecomposition::getCellContents out of bounds.\n");
      return IndexList();
    }
#endif
    return cell[i].point;
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

  int countPoints() const {
    int num = 0;
    for (int j = 0; j < size; j++) {
      num += cell[j].point.length();
    }
    return num;
  }

  IndexList neighborhood(Vector3 r) const {
    int j = getCell(r);
    int nNeighs = cell[j].neigh.length();
    IndexList ret;

    // Get the point indices in all neighbors of the home cell (inclusive).
    for (int n = 0; n < nNeighs; n++) {
      int neigh = cell[j].neigh.get(n);
      
      ret.add(cell[neigh].point);
    }
    return ret;
  }

protected:
  double cut; // cutoff length
  Matrix3 box; // periodic box of the system
  SpaceCell* cell; // contains cell attributes

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

private:  

  //CellDecomposition(const CellDecomposition&) {}
  //CellDecomposition operator=(const CellDecomposition&) {}
};


///////////////////////////////////////////////////////////////////////
// Generate a third force grid.
class ThirdForceGrid : public Grid {
public:
  ThirdForceGrid(Matrix3 box0, Vector3 origin0, double radius0, double sigma0, double dx) : Grid(box0, origin0, dx) {
    radius = radius0;
    sigma = sigma0;
    
    box = box0;
  }

  // Compute the grid using a cell decomposition, computing the energy between
  // grid points and those source in the neighborhood of the gird point.
  void compute(const Scatter& source, const CellDecomposition& decomp) {
    int j = 0;

    //IndexList one = decomp.neighborhood(Vector3(0.0));
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  Vector3 r = basis.transform(Vector3(ix, iy, iz)) + origin;
	  val[j] = localEnergy(r, source, decomp.neighborhood(r));
	  j++;
	}
      }
      if (ix % 5 == 0) printf("%2.1f percent complete\n", (100.0*j)/size);
    }
  }

  // Compute the energy at position r due to specified source points.
  double localEnergy(Vector3 r, const Scatter& source, const IndexList& neigh) const {
    const int n = neigh.length();
    //printf("%s\n\n",(const char*)neigh.toString());
   
    double e = 0.0;
    for (int i = 0; i < n; i++) {
      Vector3 r0 = source.get(neigh.get(i));
      // Transform the vector between the grid node and the source point into lattice space.
      Vector3 l = basisInv.transform(r - r0);
      // Wrap it according to the periodic boundaries.
      if (l.x < 0.5*nx) {l.x += nx;}
      if (l.x > 0.5*nx) {l.x -= nx;}
      if (l.y < 0.5*ny) {l.y += ny;}
      if (l.y > 0.5*ny) {l.y -= ny;}
      if (l.z < 0.5*nz) {l.z += nz;}
      if (l.z > 0.5*nz) {l.z -= nz;}
      
      // Transform the wrapped vector back into world space.
      Vector3 d = basis.transform(l);
      
      // Acculumulate the energy.
      e += energy(d.length());
    }

    //printf("%.10g\n", e);
    return e;
  }

  // Compute the grid by computing energy between all grid points and all sources :(.
  void compute(const Scatter& source) {
    int j = 0;

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	for (int iz = 0; iz < nz; iz++) {
	  Vector3 r = basis.transform(Vector3(ix, iy, iz)) + origin;
	  val[j] = globalEnergy(r, source);

	  j++;
	}
      }
      if (ix % 5 == 0) printf("%2.1f percent complete\n", (100.0*j)/size);
    }
  }
 
  // Compute the energy at position r due to all source points.
  // This can be slow for large systems.
  double globalEnergy(Vector3 r, const Scatter& source) const {
    const int n = source.length();
    
    double e = 0.0;
    for (int i = 0; i < n; i++) {
      Vector3 r0 = source.get(i);
      // Transform the vector between the grid node and the source point into lattice space.
      Vector3 l = basisInv.transform(r - r0);
      // Wrap it according to the periodic boundaries.
      if (l.x < 0.5*nx) {l.x += nx;}
      if (l.x > 0.5*nx) {l.x -= nx;}
      if (l.y < 0.5*ny) {l.y += ny;}
      if (l.y > 0.5*ny) {l.y -= ny;}
      if (l.z < 0.5*nz) {l.z += nz;}
      if (l.z > 0.5*nz) {l.z -= nz;}
      
      // Transform the wrapped vector back into world space.
      Vector3 d = basis.transform(l);
      
      // Acculumulate the energy.
      e += energy(d.length());
    }
    return e;
  }
   
  double energy(double r) const {
    if (r < radius) return radius - r + 0.5*sigma;
    if (r > radius + sigma) return 0.0;
    
    double dr = r - radius;
    return dr*(0.5*dr/sigma - 1.0) + 0.5*sigma;
  }

private:
  double radius, sigma;
  Matrix3 box; // periodic box of the system
  //ThirdForceGrid(const ThirdForceGrid&) {}
  //ThirdForceGrid operator=(const ThirdForceGrid&) {}
};

///////////////////////////////////////////////////////////////////////
// Histogram class.
// Inherits from Grid.
class HistogramGrid : public Grid {
public:
  HistogramGrid(const Scatter& source, int nx0, int ny0, int nz0) {
    // Get the bounds of the data.
    Vector3 r0 = source.minBound();
    Vector3 r1 = source.maxBound();
    Vector3 b = r1 - r0;

    nx = abs(nx0);
    ny = abs(ny0);
    nz = abs(nz0);
    basis = Matrix3(b.x/nx, b.y/ny, b.z/nz);
    basisInv = basis.inverse();
    origin = r0;
    size = nx*ny*nz;
    val = new double[size];
    zero();

    // Add the histogram data.
    add(source);
  }

  HistogramGrid(Matrix3 box0, Vector3 origin0, int nx, int ny, int nz) : Grid(box0, origin0, nx, ny, nz) {   
  }

  HistogramGrid(const Grid& g) : Grid(g) {   
  }

  double total() const {
    double count = 0.0;
    for (int i = 0; i < size; i++) count += val[i];
    return count;
  }

  int add(const Vector3& r) {
    int j = getCell(r);
    if (j > 0) {
      val[j] += 1.0;
      return 1;
    }
    return 0;
  }

  int add(const Scatter& pos) {
    int count = 0;
    const int n = pos.length();

    for (int i = 0; i < n; i++) count += add(pos.get(i));
    return count;
  }

  virtual Vector3 getBinCenter(int j) const {
    int iz = j%nz;
    int iy = (j/nz)%ny;
    int ix = j/(nz*ny);

    return basis.transform(Vector3(ix+0.5, iy+0.5, iz+0.5)) + origin;
  }
  
  virtual Vector3 indexToPos(int j) const {
    int iz = j%nz;
    int iy = (j/nz)%ny;
    int ix = j/(nz*ny);
    
    return Vector3(ix, iy, iz);
  }

  virtual int posToIndex(int ix, int iy, int iz) {
    return iz + iy*nz + ix*nz*ny;
  }

  // Find the cell index for a given position.
  virtual int getCell(const Vector3& r) const {
    // Transform.
    Vector3 l = basisInv.transform(r-origin);
    // Wrap into the home box.
    int cx = int(floor(l.x));
    int cy = int(floor(l.y));
    int cz = int(floor(l.z));
    
    if (cx < 0 || cx >= nx) return -1;
    if (cy < 0 || cy >= ny) return -1;
    if (cz < 0 || cz >= nz) return -1;

    return cz + cy*nz + cx*nz*ny;
  }
};

class UmbrellaWindow {
public:
  UmbrellaWindow() {
    steps = 0;
  }

  UmbrellaWindow(const Scatter& source, Vector3 center0, Vector3 spring0) {
    steps = source.length();
    center = center0;
    spring = spring0;
  }

public:
  Vector3 center;
  Vector3 spring;
  int steps;
};

class BinDim {
public:
  BinDim(Matrix3 box, Vector3 origin0, int nx0, int ny0, int nz0) {
    basis = Matrix3(box.ex()/nx0, box.ey()/ny0, box.ez()/nz0);
    origin = origin0;
    nx = nx0;
    ny = ny0;
    nz = nz0;
  }
  
  BinDim(Vector3 r0, Vector3 r1, int nx0, int ny0, int nz0) {
    origin = r0;
    Vector3 d = r1 - r0;
    basis = Matrix3(d.x/nx0, d.y/ny0, d.z/nz0);
    nx = nx0;
    ny = ny0;
    nz = nz0;
  }

   BinDim(const Scatter** win, int winN, int nx0, int ny0, int nz0) {
    Vector3 b;
    Vector3 boundMin = win[0]->minBound();
    Vector3 boundMax = win[0]->maxBound();
    for (int w = 1; w < winN; w++) {
      b = win[w]->minBound();
      if (b.x < boundMin.x) boundMin.x = b.x;
      if (b.y < boundMin.y) boundMin.y = b.y;
      if (b.z < boundMin.z) boundMin.z = b.z;

      b = win[w]->maxBound();
      if (b.x > boundMax.x) boundMax.x = b.x;
      if (b.y > boundMax.y) boundMax.y = b.y;
      if (b.z > boundMax.z) boundMax.z = b.z;
    }
    
    nx = nx0;
    ny = ny0;
    nz = nz0;
    origin = boundMin;
    Matrix3 box = Matrix3(boundMax.x-boundMin.x, boundMax.y-boundMin.y, boundMax.z-boundMin.z); 
    basis = Matrix3(box.ex()/nx0, box.ey()/ny0, box.ez()/nz0);
  }

  int length() const {
    return nx*ny*nz;
  }

public:
  int nx, ny, nz;
  Vector3 origin;
  Matrix3 basis;
};

double springEnergy(Vector3 spring, Vector3 r, Vector3 r0) {
  Vector3 d = r - r0;
  return 0.5*(spring.x*d.x*d.x + spring.y*d.y*d.y + spring.z*d.z*d.z);
}

void writeArray(const char* outFile, double* d, int n) {
  FILE* out = fopen(outFile, "w");
  for (int i = 0; i < n; i++) {
    fprintf(out, "%.10g\n", d[i]);
  } 
}

// Compute the PMF from a probability distribution.
// Write the result and a mask with the valid points.
void writePmf(const char* fileName, const char* comments, const Grid& prob) {
  int b;
  double p;
  const int n = prob.getSize();

  char maskName[256];
  sprintf(maskName, "valid_%s", fileName);

  // Find the minimum and maximum of the PMF.
  for (b = 0; b < n; b++) {
    p = prob.getValue(b);
    if (p > 0.0) break;
  }
  if (b == n) {
    printf("Error! Probability distribution is zero.\n");
    return;
  }
  double pmfMin = -log(p);
  double pmfMax = -log(p);
  for (b = 0; b < n; b++) {
    if (prob.getValue(b) > 0.0) {
      double u = -log(prob.getValue(b));
      if (u < pmfMin) pmfMin = u;
      if (u > pmfMax) pmfMax = u;
    }
  }
  
  // Calculate the PMF.
  Grid pmf(prob);
  Grid mask(prob);
  // The PMF is associated with the bin centers.
  pmf.shiftToCenters(); 
  mask.shiftToCenters();
  for (b = 0; b < n; b++) {
    if (prob.getValue(b) > 0.0) {
      double u = -log(prob.getValue(b));
      pmf.setValue(b, u - pmfMin);
      mask.setValue(b, 1);
    } else {
      pmf.setValue(b, pmfMax - pmfMin);
      mask.setValue(b, 0);
    }
  }

  pmf.write(fileName, comments);
  mask.write(maskName, comments);
}

// Omit windows with no data or large average energies.
// Deallocate their memory.
int filterWindows(Scatter** sim, UmbrellaWindow* win, int winN0, double maxWinEnergy) {
  int wg = 0;
  for (int w = 0; w < winN0; w++) {
    // Compute the average energy for each window.
    double energy = 0.0;
    for (int i = 0; i < sim[w]->length(); i++)
      energy += springEnergy(win[w].spring, win[w].center, sim[w]->get(i));
    energy = energy/sim[w]->length();

    // Omit windows that have an energy that is too large.
    if (energy < maxWinEnergy && sim[w]->length() > 0) {
      sim[wg] = sim[w];
      win[wg] = win[w];
      wg++;
      printf("Average energy for window at %s is %g kT.\n", win[w].center.toString().val(), energy);
    } else {
      printf("Average energy for window at %s is %g kT. Omitting window of %i points.\n", win[w].center.toString().val(), energy, sim[w]->length());
      delete sim[w];
    }
  }
  const int winN = wg;

  printf("%d of %d windows contribute.\n", winN, winN0);
  return winN;
}

// Make histogram for each window and put it in hist. Windows that do not contribute to the histogram will be omitted.
void histogramWindows(HistogramGrid** hist, const Scatter** sim, const UmbrellaWindow* win, int winN, BinDim bin0) {
  int w;
  // Form the bins to cover all included windows.
  BinDim bin(sim, winN, bin0.nx, bin0.ny, bin0.nz);

  // Histogram the data to form the biased probability distributions for each window.
  printf("Histogramming data...\n");
  for (w = 0; w < winN; w++) {
    hist[w] = new HistogramGrid(bin.basis, bin.origin, bin.nx, bin.ny, bin.nz);
    
    int histCount = hist[w]->add(*sim[w]);
    printf("Histogrammed %d of %d points for window %s.\n", histCount, sim[w]->length(), win[w].center.toString().val());
  }
  printf("Done.\n");
}

// Run the Weighted Histogram Analysis Method in three dimensions.
// winIndex contains the indices of the valid windows in "win".
void wham3d(const char* outName, HistogramGrid** hist, const UmbrellaWindow* win, int winN, double tol) {
  int i, w, b;

  char outFile[256];
  sprintf(outFile, "pmf_%s.dx", outName);
  char tempFile[256];
  sprintf(tempFile, "preview_%s.dx", outName);
  char histFile[256];
  sprintf(histFile, "hist_%s.dx", outName);

  
  printf("\nComputing the PMF distribution in three dimensions by the WHAM method with %d windows.\n", winN);
  printf("Based on B. Roux. Computer Physics Communications 91, 275 (1995).\n");
  
  // Get the histogram counts.
  double* histCount = new double[winN];
  for (w = 0; w < winN; w++) histCount[w] = hist[w]->total();

  // Write the total histogram.
  Grid total(*hist[0]);
  total.zero();
  for (w = 0; w < winN; w++) total.addGrid(*hist[w]);
  total.shiftToCenters();
  total.write(histFile, "total histogram");
  printf("Wrote the total histogram to %s.\n", histFile);
 
  // Scale the histograms to make them probability distributions.
  const int binN = hist[0]->length();
  const double binVol = hist[0]->getCellVolume();  
  for (w = 0; w < winN; w++) hist[w]->scale(1.0/(histCount[w]*binVol));
  
  // Compute the numerator of Roux Eq. 8 beforehand.
  printf("\nPrecomputing the numerator of Roux Eq. 8...\n");
  double* binNumer = new double[binN];
  for (b = 0; b < binN; b++) {
    double numer = 0.0;
    
    for(w = 0; w < winN; w++) {
      // Get the biased prob. distrib. at this bin center for this window.
      double p = hist[w]->getValue(b);

      // Add to the numerator sum in Roux Eq. 8.
      numer += histCount[w]*p;
    }
    binNumer[b] = numer;
  }
  printf("Done.\n");
  //writeArray("numer_C.txt", binNumer, binN);

  // Compute exp(-uBias) beforehand for Roux Eq. 8 and 9.
  printf("\nPrecomputing exp(-uBias) for Roux Eq. 8 and 9...\n");
  double** biasFactor = new double*[binN];
  for (b = 0; b < binN; b++) {
    biasFactor[b] = new double[winN];
    Vector3 binCen = hist[0]->getBinCenter(b);
    for (w = 0; w < winN; w++) {
      double uBias = springEnergy(win[w].spring, win[w].center, binCen);
      biasFactor[b][w] = exp(-uBias);
      double d = (binCen-win[w].center).length();
    }
  }
  printf("Done.\n");

  // Make the initial guess for the f constants.
  double* winF = new double[winN];
  double* winFOld = new double[winN];
  for (w = 0; w < winN; w++) {
    winF[w] = 0.0;
    winFOld[w] = 0.0;
  }

  // This contains the estimate of the unbiased prob. distrib.
  Grid prob(*hist[0]);
  prob.zero();

  //////////////////////////////
  //
  // Iterate using the WHAM equations.
  printf("\nBeginning WHAM iterations.\n");
  double maxChange = 1.0;
  int iter = 0;
  while (maxChange > tol) {
    //////////////////////////////
    // Roux Eq. 8.
    //
    // Estimate the unbiased probability distribution.
    for (b = 0; b < binN; b++) {
      double numer = binNumer[b];
      double denom = 0.0;
      // Loop through the windows to do the sum.
      for (w = 0; w < winN; w++)
	denom += histCount[w]*exp(winF[w])*biasFactor[b][w];
      prob.setValue(b, numer/denom);
    }
    
    //////////////////////////////
    // Roux Eq. 9
    //
    // Integrate to obtain the new f constants for each window.
    double* ptr = winFOld;
    winFOld = winF;
    winF = ptr;
    double sum;
    for (w = 0; w < winN; w++) {
      // Numerically integrate over the reaction coordinates (over the bins).
      sum = 0.0;
      for (b = 0; b < binN; b++)
	sum += binVol*biasFactor[b][w]*prob.getValue(b);
      winF[w] = -log(sum);
    }
 
    // Check for convergence.
    w = 1;
    double dfOld = winFOld[w] - winFOld[w-1];
    double df = winF[w] - winF[w-1];
    maxChange = fabs(df-dfOld);
    for (w = 2; w < winN; w++) {
      dfOld = winFOld[w] - winFOld[w-1];
      df = winF[w] - winF[w-1];
      if (fabs(df-dfOld) > maxChange) {
	maxChange = fabs(df-dfOld);
      }
    }
    printf("WHAM iteration %d: %.5g\n", iter, maxChange);
    iter++;

    // Write the intermediate results.
    if (iter % 100 == 0) {
      char comm[256];
      sprintf(comm, "PMF iteration %d", iter);
      writePmf(tempFile, comm, prob);
    }
  }

  // Write the F constants.
  char outNameF[256];
  sprintf(outNameF, "fconst_%s.dat", outName);
  FILE* outF = fopen(outNameF, "w");
  for (w = 0; w < winN; w++) fprintf(outF, "%.10g\n", winF[w]);
  fclose(outF);

  // We are finished with the WHAM iterations.
  printf("\nWHAM iterations finished.\n");

  // Write the PMF.
  writePmf(outFile, "PMF", prob);
  printf("Wrote %s.\n", outFile);

  // Clean up the memory.
  for (b = 0; b < binN; b++)  delete biasFactor[b];

  delete[] binNumer;
  delete[] biasFactor;
  delete[] winF;
  delete[] winFOld;
}

void readIndexFile(const char* fileName, int num, String* prefix, UmbrellaWindow* win) {
  int nRead;
  int n = 0;
  double x, y, z, kx, ky ,kz;
  char line[256];
  char ind[256];

  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("readIndexFile Couldn't open file %s\n.",fileName);
    exit(-1);
  }

  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead=sscanf(line,"%s %lf %lf %lf %lf %lf %lf",ind,&x,&y,&z,&kx,&ky,&kz);
    if (nRead >= 7) {
      prefix[n].add(ind);
      win[n].center.x = x;
      win[n].center.y = y;
      win[n].center.z = z;
      win[n].spring.x = kx;
      win[n].spring.y = ky;
      win[n].spring.z = kz;
      n++;
      if (n >= num) break;
    } else {
      printf("Warning! Improperly formatted index file %s.\n", fileName);
      printf("Simulation index file format: winDataFile x0 y0 z0 kx ky kz.\n");
    }
  }
    
  fclose(inp);
}

int countIndexFile(const char* fileName) {
  int nRead;
  int n = 0;
  double x, y, z, kx, ky, kz;
  char line[256];
  int ind;

  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("countIndexFile Couldn't open file %s\n.",fileName);
    exit(-1);
  }

  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead=sscanf(line,"%s %lf %lf %lf %lf %lf %lf",&ind,&x,&y,&z,&kx,&ky,&kz);
    if (nRead >= 7) n++;
  }
    
  fclose(inp);
  return n;
}

BinDim readBinDimensions(const char* fileName) {
  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("readBinDimensions Couldn't open file %s\n.",fileName);
    exit(-1);
  }
  int nRead;
  double x0, y0, z0, x1, y1, z1;
  int nx, ny, nz;

  char line[256];
  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead=sscanf(line, "%lf %lf %lf %lf %lf %lf %d %d %d",&x0,&y0,&z0,&x1,&y1,&z1,&nx,&ny,&nz);
    if (nRead != 9) {
      printf("readBinDimensions File %s improperly formatted.\n", fileName);
      printf("Format is 'x0 y0 z0 x1 y1 z1 nx ny nz'.\n");
      exit(1);
    }
  }

  return BinDim(Vector3(x0,y0,z0), Vector3(x1,y1,z1), nx, ny, nz);
}

///////////////////////////////////////////////////////////////////////
// Drivers
// indexFile has records "prefix x y z kx ky kz".
// Spring constants are assumed to be in kcal/mol A^-2.
int main(int argc, char* argv[]) {
  if ( argc != 6 ) {
    printf("Usage: %s indexFile binDimensionsFile cutTime0 cutTime1 outName\n", argv[0]);
    return 0;
  }

  // Parameters:
  double tol = 0.001;
  //Vector3 springK(0.8, 0.8, 4.0); // in kcal/mol A^-2
  const double cutTime0 = strtod(argv[3], NULL);
  const double cutTime1 = strtod(argv[4], NULL);
  // Input:
  const char* indexFile = argv[1];
  const char* binDimFile = argv[2];
  // Output:
  const char* outName = argv[5];
  
  const double kBT = 0.5862292; // in kcal/mol at 295 K
  // Windows with average energies greater than this are discarded.
  const double maxWinEnergy = 3.0*0.8; // in kT
  int w;

  // Read the bin dimensions from a file.
  //BinDim bin(win, winN, nx, ny, nz);
  BinDim bin = readBinDimensions(binDimFile);
  printf("Read bin dimensions from %s.\n", binDimFile);

  int winN = countIndexFile(indexFile);
  printf("Found %d windows in simulation index file %s.\n", winN, indexFile);
  
  // Window variables
  String* winDataFile = new String[winN];
  UmbrellaWindow* win = new UmbrellaWindow[winN];
  readIndexFile(indexFile, winN, winDataFile, win);

  // Convert to k_B T/A^2.
  for (w = 0; w < winN; w++) win[w].spring /= kBT;
  printf("Converted the spring constants from kcal/mol/x^2 to k_B T/x^2.\n");

  // Read the simulation data.
  printf("Reading simulation data, ignoring data with %g < time < %g.\n", cutTime0, cutTime1);
  Scatter** sim = new Scatter*[winN];
  char dataFile[256];
  for (w = 0; w < winN; w++) {
    sim[w] = new Scatter(winDataFile[w], cutTime0, cutTime1);
    win[w].steps = sim[w]->length();

    printf("Read %d data points from %s.\n", win[w].steps, winDataFile[w].val());
  }

  // Filter the windows.
  const int winN1 = filterWindows(sim, win, winN, maxWinEnergy);

  // Make the histograms.
  HistogramGrid** hist = new HistogramGrid*[winN1];
  histogramWindows(hist, (const Scatter**)sim, win, winN1, bin);

  // Run WHAM.
  wham3d(outName, hist, win, winN1, tol);

  // Clean up.
  for (w = 0; w < winN1; w++) delete sim[w];
  for (w = 0; w < winN1; w++)  delete hist[w];
 
  delete[] sim;
  delete[] win;
  delete[] winDataFile;
  delete[] hist;

  return 0;
}
