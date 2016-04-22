function [data grid delta origin] = readDxData(fileName)
% Read data formatted as:
% nx ny nz ox oy oz dxx dyx dzx dxy dyy dzy dxz dyz dzz val0 val1 val2 ...
% [DATA GRID DELTA ORIGIN] = READDX(FILENAME) returns the potential
% value at each point in DATA, the number of grid points along each
% cell axis in GRID (3-vector), the lattice vectors in DELTA
% (3x3 matrix), and the cell origin in ORIGIN (3-vector).
% Data order is z fast, y medium, and x slow.
% Author: Jeff Comer <jcomer2@illinois.edu>

in = fopen(fileName, 'r');
if in == -1
    error('readDxData:fileNotFound', 'The .dx file was not found.')
end
disp(sprintf('Reading file %s.', fileName))

% Read the file contents.
stuff = dlmread(fileName, ' ');
grid = stuff(1:3);
origin = stuff(4:6);
delta = [stuff(7:9) stuff(10:12) stuff(13:15)];
data = stuff(16:end);
disp(sprintf('Successfully read %i values.', length(data)))





