
///////////////////////////////////////////////////////////////////////
// Cell decomposition of points.
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef HISTOGRAMGRID_H
#define HISTOGRAMGRID_H

#include "Grid.H"

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

#endif
