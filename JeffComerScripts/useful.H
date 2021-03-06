///////////////////////////////////////////////////////////////////////
// Useful classes.
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef USEFUL_H
#define USEFUL_H

#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
using namespace std;

bool isReal(char c) {
  char num[] = "0123456789-eE.";
  
  for (int i = 0; i < 14; i++) {
    if (c == num[i]) return true;
  }
  return false;
}

bool isInt(char c) {
  char num[] = "0123456789-";
  
  for (int i = 0; i < 11; i++) {
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
    len = s.len;
    cap = s.len;
    c = new char[cap];
    for (int i = 0; i < s.len; i++) c[i] = s.c[i];
    return *this;
  }
  
  void add(char c0) {
    int n = 2;
    len--;

    if (n + len > cap) grow(n + len);
    c[len] = c0;
    c[len+1] = '\0';
    len += n;
  }

  void add(const char* s) {
    const int n = strlen(s) + 1;
    len--;

    if (n + len > cap) grow(n + len);
    for (int i = 0; i < n; i++) c[i+len] = s[i];
    len += n;
  }

  void add(String& s) {
    len--;

    if (len + s.len > cap) grow(len + s.len);
    for (int i = 0; i < s.len; i++) c[i+len] = s.c[i];
    len += s.len;
  }
  
  int length() const {
    return len-1;
  }

  // Negative indices go from the end.
  String range(int first, int last) const {
    String ret;
    const int len1 = len - 1;

    first = first % len1;
    last = last % len1;

    if (first < 0) first += len1;
    if (last < 0) last += len1;
    for(int i = first; i <= last; i++) ret.add(c[i]);

    return ret;
  }

  String trim() {
    int i0, i1;
    for (i0 = 0; i0 < len; i0++) if (!isWhite(c[i0])) break;
    for (i1 = len-1; i1 >= 0; i1++) if (!isWhite(c[i0])) break;
    return range(i0,i1);
  }
  
  static bool isWhite(char c) {
    char num[] = " \n\t\v\b\r\f\a";

    for (int i = 0; i < 8; i++) {
      if (c == num[i]) return true;
    }
    return false;
  }

  int tokenCount() const {
    const int len1 = len - 1;

    int count = 0;

    // Skip leading whitespace.
    int i = 0;
    while (isWhite(c[i]) && i < len1) i++;

    if (i >= len1) return 0;
    
    // Find the tokens.
    bool onToken = true;
    count++;
    while (i < len1) {
      if (onToken) {
	if (isWhite(c[i])) onToken = false;
      } else {
	if (!isWhite(c[i])) {
	  onToken = true;
	  count++;
	}
      }
      i++;
    }

    return count;
  }

  int tokenize(String* tokenList) const {
    const int len1 = len - 1;
    const int n = tokenCount();

    // Skip leading whitespace.
    int i = 0;
    while (isWhite(c[i]) && i < len1) i++;

    // Find the tokens.
    bool onToken = true;
    int j = 0;
    int tokenStart = i;
    while (i < len1) {
      if (onToken) {
	if (isWhite(c[i])) {
	  onToken = false;
	  tokenList[j] = range(tokenStart,i-1);
	  //printf("(%s)\n", tokenList[j].val());
	  j++;
	}
      } else {
	if (!isWhite(c[i])) {
	  onToken = true;
	  tokenStart = i;
	}
      }
      i++;
    }

    // Grab the last token if there is one.
    if (onToken) tokenList[n-1] = range(tokenStart,-1);
    return n;
  }

  int tokenCount(char delim) const {
    const int len1 = len - 1;

    int count = 0;

    // Skip leading whitespace.
    int i = 0;
    while (c[i] == delim && i < len1) i++;

    if (i >= len1) return 0;
    
    // Find the tokens.
    bool onToken = true;
    count++;
    while (i < len1) {
      if (onToken) {
	if (c[i] == delim) onToken = false;
      } else {
	if (c[i] != delim) {
	  onToken = true;
	  count++;
	}
      }
      i++;
    }

    return count;
  }

  int tokenize(String* tokenList, char delim) const {
    const int len1 = len - 1;
    const int n = tokenCount(delim);

    // Skip leading whitespace.
    int i = 0;
    while (c[i] == delim && i < len1) i++;

    // Find the tokens.
    bool onToken = true;
    int j = 0;
    int tokenStart = i;
    while (i < len1) {
      if (onToken) {
	if (c[i] == delim) {
	  onToken = false;
	  tokenList[j] = range(tokenStart,i-1);
	  //printf("(%s)\n", tokenList[j].val());
	  j++;
	}
      } else {
	if (c[i] != delim) {
	  onToken = true;
	  tokenStart = i;
	}
      }
      i++;
    }

    // Grab the last token if there is one.
    if (onToken) tokenList[n-1] = range(tokenStart,-1);
    return n;
  }

  ~String() {
    delete[] c;
  }

  void lower() {
    for(int i = 0; i < len; i++) {
      if (c[i] >= 0x41 && c[i] <= 0x5A) c[i] += 0x20;
    }
  }
  
  void upper() {
    for(int i = 0; i < len; i++) {
      if (c[i] >= 0x61 && c[i] <= 0x7A) c[i] -= 0x20;
    }
  }

  bool operator==(const String& s) const {
    if (len != s.len) return false;
    for (int i = 0; i < len; i++) if (c[i] != s.c[i]) return false;
    return true;
  }

  bool operator!=(const String& s) const {
    return !(*this == s);
  }

  operator const char*() const {
    return c;
  }

  const char* val() const {
    return c;
  }

  String getNumbers() const {
    String ret;
    for(int i = 0; i < len; i++) {
      if (isInt(c[i])) ret.add(c[i]);
    }
    return ret;
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

  Vector3(double s):x(s), y(s), z(s) {
  }

  Vector3(const Vector3& v):x(v.x), y(v.y), z(v.z) {
  }
  
  Vector3(double x0, double y0, double z0) {
    x = x0;
    y = y0;
    z = z0;
  }

  Vector3(const double* d) {
    x = d[0];
    y = d[1];
    z = d[3];
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
    v.x = y*w.z - z*w.y;
    v.y = z*w.x - x*w.z;
    v.z = x*w.y - y*w.x;
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

  Matrix3(const double* d) {
    exx = d[0];
    exy = d[1];
    exz = d[2];
    eyx = d[3];
    eyy = d[4];
    eyz = d[5];
    ezx = d[6];
    ezy = d[7];
    ezz = d[8];
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


  const Matrix3 operator*(const Matrix3& m) const {
    Matrix3 ret;
    ret.exx = exx*m.exx + exy*m.eyx + exz*m.ezx;
    ret.eyx = eyx*m.exx + eyy*m.eyx + eyz*m.ezx;
    ret.ezx = ezx*m.exx + ezy*m.eyx + ezz*m.ezx;

    ret.exy = exx*m.exy + exy*m.eyy + exz*m.ezy;
    ret.eyy = eyx*m.exy + eyy*m.eyy + eyz*m.ezy;
    ret.ezy = ezx*m.exy + ezy*m.eyy + ezz*m.ezy;

    ret.exz = exx*m.exz + exy*m.eyz + exz*m.ezz;
    ret.eyz = eyx*m.exz + eyy*m.eyz + eyz*m.ezz;
    ret.ezz = ezx*m.exz + ezy*m.eyz + ezz*m.ezz;
    return ret;
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

  
  Matrix3 transform(const Matrix3& m) const {
    Matrix3 ret;
    ret.exx = exx*m.exx + exy*m.eyx + exz*m.ezx;
    ret.eyx = eyx*m.exx + eyy*m.eyx + eyz*m.ezx;
    ret.ezx = ezx*m.exx + ezy*m.eyx + ezz*m.ezx;

    ret.exy = exx*m.exy + exy*m.eyy + exz*m.ezy;
    ret.eyy = eyx*m.exy + eyy*m.eyy + eyz*m.ezy;
    ret.ezy = ezx*m.exy + ezy*m.eyy + ezz*m.ezy;

    ret.exz = exx*m.exz + exy*m.eyz + exz*m.ezz;
    ret.eyz = eyx*m.exz + eyy*m.eyz + eyz*m.ezz;
    ret.ezz = ezx*m.exz + ezy*m.eyz + ezz*m.ezz;
    return ret;
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

  int find(int key) {
    for(int i = 0; i < num; i++) if (lis[i] == key) return i;
    return -1;
  }

private:
  int num;
  int maxnum;
  int* lis;
};
#endif
