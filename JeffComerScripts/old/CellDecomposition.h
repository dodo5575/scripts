///////////////////////////////////////////////////////////////////////
// Cell decomposition of points.
// Author: Jeff Comer <jcomer2@illinois.edu>
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
