///////////////////////////////////////////////////////////////////////
// Grid base class.
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef GRID_H
#define GRID_H

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

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
    
      if (isInt(line[0]) && n < size) {
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
	for (int i = 0; i < 5; i++) start[i] = line[i];
	start[5] = '\0';

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
	  //printf("%s", line);
	  // Get the system dimensions.
	  if (line[7] != '1') continue;
	  int read = sscanf(line, "object 1 class gridpositions counts %d %d %d\n", &nx, &ny, &nz);
	  //printf("Size: %d %d %d\n", nx, ny, nz);
	  if (read == 3) {
	    size = nx*ny*nz;
	    val = new double[size];
	    zero();
	  }
	}
      }
    }
    fclose(inp);

    basis = Matrix3(base[0], base[1], base[2]);
    basisInv = basis.inverse();
    if (size == 0) {
      printf("Error Grid:Grid Improperly formatted dx file %s.\n",fileName);
      exit(-1);
    }
  }

  // Writes the grid as a file in the dx format.
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
      fprintf(out, "%.10g %.10g %.10g\n", val[i], val[i+1], val[i+2]);
    }
    if (mod == 1) {
      fprintf(out, "%.10g\n", val[size-1]);
    } else if (mod == 2) {
      fprintf(out, "%.10g %.10g\n", val[size-2], val[size-1]);
    }
    fclose(out);
  }

  // Writes the grid data as a single column in the order:
  // nx ny nz ox oy oz dxx dyx dzx dxy dyy dzy dxz dyz dzz val0 val1 val2 ...
  virtual void writeData(const char* fileName) {
    // Open the file.
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s\n.",fileName);
      exit(-1);
    }

    fprintf(out, "%.10g\n%.10g\n%.10g\n", nx, ny, nz);
    fprintf(out, "%.10g\n%.10g\n%.10g\n", origin.x, origin.y, origin.z);
    fprintf(out, "%.10g\n%.10g\n%.10g\n", basis.exx, basis.eyx, basis.ezx);
    fprintf(out, "%.10g\n%.10g\n%.10g\n", basis.exx, basis.eyx, basis.ezx);
    fprintf(out, "%.10g\n%.10g\n%.10g\n", basis.exx, basis.eyx, basis.ezx);

    for (int i = 0; i < size; i++) fprintf(out, "%.10g\n", val[i]);
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
  virtual Vector3 getDestination() const {
    return basis.transform(Vector3(nx,ny,nz)) + origin;
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
      fprintf(out, "%0.10g %0.10g\n", z, v);
    }
    fclose(out);
  }

  // Get the average value of a section of the grid.
  virtual double averageSection(int iz) {
    const int nynz = ny*nz;
    double sum = 0.0;

    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	int j = iz + iy*nz + ix*nynz;
	sum += val[j];
      }
    }
      
    double v = sum/(nx*ny);
    return v;
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
      fprintf(out, "%0.10g %0.10g\n", z, v);
    }

    fclose(out);
  }


  // Compute the average profile along z.
  virtual void averageProfileZBoltzmann(const char* fileName) {
    const int nynz = ny*nz;
    
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s\n.",fileName);
      exit(-1);
    }

   
    for (int iz = 0; iz < nz; iz++) {
      double sumNum = 0.0;
      double sumDen = 0.0;

      for (int ix = 0; ix < nx; ix++) {
	for (int iy = 0; iy < ny; iy++) {
	  int j = iz + iy*nz + ix*nynz;
	  double weight = exp(-val[j]);

	  sumNum += val[j]*weight;
	  sumDen + weight;
	}
      }
      
      double v = sumNum/sumDen;
      double z = origin.z + iz*basis.ezz;
      fprintf(out, "%0.10g %0.10g\n", z, v);
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
      fprintf(out, "%0.10g %0.10g\n", z, v);
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
    
    char profileName[256];
    sprintf(profileName, "%s.profile", fileName);
    FILE* outProfile = fopen(profileName, "w");
    if (outProfile == NULL) {
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
	  fprintf(outProfile, "%.10g ", val[i]);

	  if (val[i] < minVal) {
	    minVal = val[i];
	    minIz = iz;
	  }
	  i++;
	}

	// Print the map of the minimum values.
	Vector3 r = basis.transform(Vector3(ix,iy,minIz)) + origin;
	fprintf(out, "%.10g %.10g %.10g %.10g\n", r.x, r.y, r.z, minVal);
	fprintf(outProfile, "\n");
      }
    }
    fclose(outProfile);
    fclose(out);
  }

  double interpolatePotential(Vector3 pos) {
    // Find the home node.
    Vector3 l = basisInv.transform(pos - origin);
    int homeX = int(floor(l.x));
    int homeY = int(floor(l.y));
    int homeZ = int(floor(l.z));
    
    // Get the array jumps with shifted indices.
    int jump[3];
    jump[0] = nz*ny;
    jump[1] = nz;
    jump[2] = 1;

    // Shift the indices in the home array.
    int home[3];
    home[0] = homeX;
    home[1] = homeY;
    home[2] = homeZ;

    // Shift the indices in the grid dimensions.
    int g[3];
    g[0] = nx;
    g[1] = ny;
    g[2] = nz;

    // Get the interpolation coordinates.
    double w[3];
    w[0] = l.x - homeX;
    w[1] = l.y - homeY;
    w[2] = l.z - homeZ;

    // Find the values at the neighbors.
    double g1[4][4][4];
    int ix, iy, iz;
    for (int ix = 0; ix < 4; ix++) {
      for (int iy = 0; iy < 4; iy++) {
	for (int iz = 0; iz < 4; iz++) {
	  // Wrap around the periodic boundaries. 
	  int jx = ix-1 + home[0];
	  jx = wrap(jx, g[0]);
	  int jy = iy-1 + home[1];
	  jy = wrap(jy, g[1]);
	  int jz = iz-1 + home[2];
	  jz = wrap(jz, g[2]);
	  
	  int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
	  g1[ix][iy][iz] = val[ind];
	}
      }
    }
    double a0, a1, a2, a3;
  
    // Mix along x, taking the derivative.
    double g2[4][4];
    for (int iy = 0; iy < 4; iy++) {
      for (int iz = 0; iz < 4; iz++) {

	a3 = 0.5*(-g1[0][iy][iz] + 3*g1[1][iy][iz] - 3*g1[2][iy][iz] + g1[3][iy][iz]);
	a2 = 0.5*(2*g1[0][iy][iz] - 5*g1[1][iy][iz] + 4*g1[2][iy][iz] - g1[3][iy][iz]);
	a1 = 0.5*(-g1[0][iy][iz] + g1[2][iy][iz]);
	a0 = g1[1][iy][iz];

	g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
      }
    }

    // Mix along y.
    double g3[4];
    for (int iz = 0; iz < 4; iz++) {
      a3 = 0.5*(-g2[0][iz] + 3*g2[1][iz] - 3*g2[2][iz] + g2[3][iz]);
      a2 = 0.5*(2*g2[0][iz] - 5*g2[1][iz] + 4*g2[2][iz] - g2[3][iz]);
      a1 = 0.5*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
    }

    // Mix along z.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    return a3*w[2]*w[2]*w[2] + a2*w[2]*w[2] + a1*w[2] + a0;
  }

  Vector3 interpolateForce(Vector3 pos) const {
    Vector3 f;
    f.x = interpolateGrad(pos, 0);
    f.y = interpolateGrad(pos, 1);
    f.z = interpolateGrad(pos, 2);
    return f;
  }

  static int wrap(int i, int n) {
    if (i < 0) {
      i %= n;
      i += n;
    } else if (i >= n) {
      i %= n;
    }
    return i;
  }

  double interpolateGrad(Vector3 pos, int dir) const {
    // Find the home node.
    Vector3 l = basisInv.transform(pos - origin);
    int homeX = int(floor(l.x));
    int homeY = int(floor(l.y));
    int homeZ = int(floor(l.z));

    // Find the indices of the other directions.
    int dirX = dir % 3;
    int dirY = (dir+1) % 3;
    int dirZ = (dir+2) % 3;

    // Get the array jumps with shifted indices.
    int jump[3];
    jump[dirX] = nz*ny;
    jump[dirY] = nz;
    jump[dirZ] = 1;
    //printf("dir %d %d %d\n", dirX, dirY, dirZ);

    // Shift the indices in the home array.
    int home[3];
    home[dirX] = homeX;
    home[dirY] = homeY;
    home[dirZ] = homeZ;

    // Shift the indices in the grid dimensions.
    int g[3];
    g[dirX] = nx;
    g[dirY] = ny;
    g[dirZ] = nz;

    // Get the interpolation coordinates.
    double w0[3];
    w0[0] = l.x - homeX;
    w0[1] = l.y - homeY;
    w0[2] = l.z - homeZ;
    double w[3];
    w[dirX] = w0[0];
    w[dirY] = w0[1];
    w[dirZ] = w0[2];

    // Find the values at the neighbors.
    double g1[4][4][4];
    int ix, iy, iz;
    for (ix = 0; ix < 4; ix++) {
      for (iy = 0; iy < 4; iy++) {
	for (iz = 0; iz < 4; iz++) {
	  // Wrap around the periodic boundaries.
	  int jx = ix-1 + home[0];
	  jx = wrap(jx, g[0]);
	  int jy = iy-1 + home[1];
	  jy = wrap(jy, g[1]);
	  int jz = iz-1 + home[2];
	  jz = wrap(jz, g[2]);
	  
	  int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
	  g1[ix][iy][iz] = val[ind];
	}
      }
    }
    double a0, a1, a2, a3;

    // Mix along x, taking the derivative.
    double g2[4][4];
    for (iy = 0; iy < 4; iy++) {
      for (iz = 0; iz < 4; iz++) {
	a3 = 0.5*(-g1[0][iy][iz] + 3*g1[1][iy][iz] - 3*g1[2][iy][iz] + g1[3][iy][iz]);
	a2 = 0.5*(2*g1[0][iy][iz] - 5*g1[1][iy][iz] + 4*g1[2][iy][iz] - g1[3][iy][iz]);
	a1 = 0.5*(-g1[0][iy][iz] + g1[2][iy][iz]);
	a0 = g1[1][iy][iz];

	g2[iy][iz] = 3.0*a3*w[0]*w[0] + 2.0*a2*w[0] + a1;
	//g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
      }
    }

    // Mix along y.
    double g3[4];
    for (iz = 0; iz < 4; iz++) {
      a3 = 0.5*(-g2[0][iz] + 3*g2[1][iz] - 3*g2[2][iz] + g2[3][iz]);
      a2 = 0.5*(2*g2[0][iz] - 5*g2[1][iz] + 4*g2[2][iz] - g2[3][iz]);
      a1 = 0.5*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
    }

    // Mix along z.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    return -(a3*w[2]*w[2]*w[2] + a2*w[2]*w[2] + a1*w[2] + a0);
  } 

protected:
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
#endif
