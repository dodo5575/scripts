///////////////////////////////////////////////////////////////////////
// Generate a third force grid.
// Author: Jeff Comer <jcomer2@illinois.edu>
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
// Driver
int mainThirdForce(int argc, char* argv[])
{
  char s[128];

  if ( argc != 7 ) {
    printf(
      "Usage: %s coordFile basisFile dx radius0 sigma outFile\n",
      argv[0]);
    return 0;
  }

  // Extract the parameters.
  double dx = strtod(argv[3],NULL);
  double radius = strtod(argv[4],NULL);
  double sigma = strtod(argv[5],NULL);
  printf("\n********************************\nThird force grid: \n");
  printf("The resolution is %.10g.\n", dx);
  printf("Creating a grid with radius0 = %.10g and sigma = %.10g...\n",
	 radius, sigma);

  // Get the system dimensions.
  Scatter sysVec(argv[2]);
  Matrix3 sys = sysVec.topMatrix();
  printf("System dimensions:\n%s\n", (const char*)sys.toString());
  Vector3 origin = (sys.ex() + sys.ey() + sys.ez())*-0.5;
  printf("System origin:%s\n", (const char*)origin.toString());

  // Get the source points.
  Scatter source(argv[1]);
  printf("\nLoaded %i source points from %s\n", source.length(), argv[1]);

  // Generate the cell decomposition.
  double cutoff = 1.8*(radius+sigma);
  printf("\nGenerating cell decomposition for %s with cutoff = %.10g...\n", argv[1], cutoff);
  CellDecomposition cell(sys, origin, cutoff);
  printf("Cell decomposition of %i cells\n", cell.length());
  cell.decompose(source);
  printf("Cell decomposition contains %i points\n", cell.countPoints());

  // Validate the cell decomposition.
  const int nCells = cell.length();
  const int nPoints = source.length();
  
  // Organize the file comments.
  char comments[128];
  //sprintf(comments, "ThirdForce grid for %s with dx = %.10g, radius0 = %.10g, and sigma = %.10g", argv[1], dx, radius, sigma);
  comments[0] = '\0';
  
  // Generate the grid.
  ThirdForceGrid third(sys, origin, radius, sigma, dx);
  printf("\nComputing the grid...\n");
  third.compute(source, cell);
  //third.compute(source);
  printf("Writing the grid...\n");
  third.write(argv[6], comments);
  printf("Wrote %s successfully.\n", argv[6]);

  return 0;
}
