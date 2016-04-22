// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "useful.H"
#include "Grid.H"

using namespace std;

// Do the first n characters match?
bool match(const char* a, const char* b, int n) {
  for (int i = 0; i < n; i++) {
    if (a[i] != b[i]) return false;
  }
  return true;
}

bool drawAll(Grid& g, const char* line) {
  double v;
  char kind[128];
  int nRead = sscanf(line, "%s %lf", kind, &v);
  if (nRead != 2) return false;

  const int n = g.length();
  for (int i = 0; i < n; i++) g.setValue(i, v);

  return true;
}

bool drawBlur(Grid& g, const char* line, bool wrapping) {
  // WRAPPING is currently ignored!

  int count;
  char kind[128];
  int nRead = sscanf(line, "%s %d", kind, &count);

  if (nRead != 2) return false;
  for (int j = 0; j < count; j++) g.blur();

  return true;
}

bool drawSphere(Grid& g, const char* line, bool wrapping) {
  double x, y, z, rad, v;
  char kind[128];
  int nRead = sscanf(line, "%s %lf %lf %lf %lf %lf", kind, &x, &y, &z, &rad, &v);

  if (nRead != 6) return false;
  Vector3 cen(x, y, z);
  double radSq = rad*rad;

  const int n = g.length();
  for (int i = 0; i < n; i++) {
    Vector3 r = g.getPosition(i);
    Vector3 d = r - cen;
    
    if (wrapping) d = g.wrapDiff(d);
    if (d.length2() <= radSq) g.setValue(i, v);
  }

  return true;
}

bool inBox(Vector3 r, Vector3 a, Vector3 b) {
  if (r.x >= a.x && r.x < b.x && r.y >= a.y && r.y < b.y && r.z >= a.z && r.z < b.z) return true;
  return false;
}

bool drawBox(Grid& g, const char* line, bool wrapping) {
  double ax, ay, az, bx, by, bz, v;
  char kind[128];
  int nRead = sscanf(line, "%s %lf %lf %lf %lf %lf %lf %lf", kind, &ax, &ay, &az, &bx, &by, &bz, &v);
  const int n = g.length();

  // Validate the box.
  if (nRead != 8) return false;
  Vector3 a(0.0);
  Vector3 b(0.0);
  if (ax < bx) { a.x = ax; b.x = bx; }
  else { a.x = bx; b.x = ax; }
  if (ay < by) { a.y = ay; b.y = by; }
  else { a.y = by; b.y = ay; }
  if (az < bz) { a.z = az; b.z = bz; }
  else { a.z = bz; b.z = az; }
    
  if (wrapping) {
    // Vectors for images.
    Vector3 image[27];
    int j = 0;
    Matrix3 box = g.getBox();
    Vector3 ex(box.ex());
    Vector3 ey(box.ey());
    Vector3 ez(box.ez());
    for (int ix = -1; ix < 1; ix++) {
      for (int iy = -1; iy < 1; iy++) {
	for (int iz = -1; iz < 1; iz++) {
	  image[j] = ex*ix + ey*iy + ez*iz;
	  j++;
	}
      }      
    }

    // Loop through all nodes.
    for (int i = 0; i < n; i++) {
      Vector3 r = g.getPosition(i);
       
      // Check all images.
      for (int j = 0; j < 27; j++) {
	Vector3 r0 = r + image[j];
	if (inBox(r0, a, b)) {
	  g.setValue(i, v);
	  break;
	}
      }
    }

  } else {
    for (int i = 0; i < n; i++) {
      Vector3 r = g.getPosition(i);
      if (inBox(r, a, b)) g.setValue(i, v);
    }
  }

  return true;
}

// Distance from a line to a point.
// u is a unit vector.
double lineDistance(Vector3 r0, Vector3 u, Vector3 r) {
  Vector3 d = r - r0;
  double comp = d.dot(u);
  return sqrt(d.length2() - comp*comp);
}

bool drawCylinder(Grid& g, const char* line, bool wrapping) {
  double x, y, z, ax, ay, az, rad, v;
  char kind[128];
  int nRead = sscanf(line, "%s %lf %lf %lf %lf %lf %lf %lf %lf", kind, &x, &y, &z, &ax, &ay, &az, &rad, &v);

  if (nRead != 9) return false;
  Vector3 cen(x, y, z);
  Vector3 axis(ax, ay, az);
  axis /= axis.length();

  const int n = g.length();
  if (wrapping) {
    cen = g.wrap(cen);

    // Vectors for images.
    Vector3 image[27];
    int j = 0;
    Matrix3 box = g.getBox();
    Vector3 ex(box.ex());
    Vector3 ey(box.ey());
    Vector3 ez(box.ez());
    for (int ix = -1; ix < 1; ix++) {
      for (int iy = -1; iy < 1; iy++) {
	for (int iz = -1; iz < 1; iz++) {
	  image[j] = ex*ix + ey*iy + ez*iz;
	  j++;
	}
      }      
    }

    // Loop through all nodes.
    for (int i = 0; i < n; i++) {
      Vector3 r = g.getPosition(i);
       
      // Check all images.
      for (int j = 0; j < 27; j++) {
	Vector3 r0 = cen + image[j];
	if (lineDistance(r0, axis, r) <= rad) {
	  g.setValue(i, v);
	  break;
	}
      }
    }

  } else {
    for (int i = 0; i < n; i++) {
      Vector3 r = g.getPosition(i);
      if (lineDistance(cen, axis, r) <= rad) g.setValue(i, v);
    }
  }

  return true;
}

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc < 4 ) {
    printf("Set regions of the grid to given values.\n");
    printf("Usage: %s srcGrid geometryFile outFile\n\n", argv[0]);
    printf("The geometryFile contains any number of commands like the following:\n");
    printf("box ax ay az bx by bz value\n");
    printf("cylinder rx ry rz normal.x normal.y normal.z radius value\n");
    printf("sphere rx ry rz radius value\n");
    printf("blur count\n");
    printf("all value\n");

    printf("\nWrapping around periodic boundaries can be switched on and off with: \n");
    printf("wrap on\n");
    printf("wrap off\n");
    printf("Wrapping is on by default.\n");

    return 0;
  }

  const char* inFile = argv[1];
  const char* geometryFile = argv[2];
  const char* outFile = argv[argc-1];
  bool wrapping = true;

  // Load the first grid.
  Grid src(inFile);
  printf("Loaded `%s'.\n", argv[1]);

  // Open the geometry file.
  FILE* inp = fopen(geometryFile,"r");
  if (inp == NULL) {
    printf("Couldn't open geometry file %s.\n", geometryFile);
    exit(-1);
  }

  char line[256];  
  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Follow commands.
    if (match("wrap on", line, 7)) {
      printf("wrap on\n");
      wrapping = true;
    } else if (match("wrap off", line, 8)) {
      printf("wrap off\n");
      wrapping = false;
    } else if (match("sphere", line, 6)) {
      printf("sphere\n");
      if (!drawSphere(src, line, wrapping)) printf("ERROR! Invalid command: %s\n", line);
    } else if (match("cylinder", line, 8)) {
      printf("cylinder\n");
      if (!drawCylinder(src, line, wrapping)) printf("ERROR! Invalid command: %s\n", line);
    } else if (match("all", line, 3)) {
      printf("all\n");
      if (!drawAll(src, line)) printf("ERROR! Invalid command: %s\n", line);
    }else if (match("box", line, 3)) {
      printf("box\n");
      if (!drawBox(src, line, wrapping)) printf("ERROR! Invalid command: %s\n", line);
    } else if (match("blur", line, 4)) {
      printf("blur\n");
      if (!drawBlur(src, line, wrapping)) printf("ERROR! Invalid command: %s\n", line);
    }
  }

  char comments[256];
  sprintf(comments, "%s modified by %s", inFile, geometryFile);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
