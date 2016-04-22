function n = cellLookup(pos, d, origin, grid)
% Organize a set of 3D points into cells.
% N = cellLookup(POS, D, ORIGIN, GRID).
% Find the cell index of the point for the cell decomposition given by
% d, origin, and grid.
% Bags go x slow, y medium, z fast.
% Author: Jeff Comer <jcomer2@illinois.edu>

ix = floor((pos(1) - origin(1))/d);
iy = floor((pos(2) - origin(2))/d);
iz = floor((pos(3) - origin(3))/d);
n = 1 + iz + iy*grid(3) + ix*grid(2)*grid(3);



