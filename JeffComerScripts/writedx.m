function writedx(data, grid, delta, origin, fileName)
% Write a .dx file.
% writedx(DATA,GRID,DELTA,ORIGIN,FILENAME) writes the
% value at each point in DATA, the number of grid points along each
% cell axis in GRID (3-vector), the lattice vectors in DELTA
% (3x3 matrix), and the cell origin in ORIGIN (3-vector) to a .dx file.
% Data order is z fast, y medium, and x slow.
% Author: Jeff Comer <jcomer2@illinois.edu>

out = fopen(fileName, 'w');
if out == -1
    error('writedx:ioError', 'Could not write the dx file.')
end
disp(sprintf('Writing file %s.', fileName))

% Write the headers.
n = length(data);
fprintf(out, 'object 1 class gridpositions counts %d %d %d\n',...
    grid(1), grid(2), grid(3));
fprintf(out, 'origin %.6f %.6f %.6f\n', origin(1), origin(2), origin(3));
for j=1:3
    fprintf(out, 'delta %.6f %.6f %.6f\n',...
        delta(1,j), delta(2,j), delta(3,j));
end
fprintf(out, 'object 2 class gridconnections counts %d %d %d\n',...
    grid(1), grid(2), grid(3));
fprintf(out,...
    'object 3 class array type double rank 0 items %d data follows\n',n);

% Write the data.
n3 = n/3;
k = 1;
for j=1:n3
    fprintf(out, '%.6f %.6f %.6f\n', data(k), data(k+1), data(k+2));
    k = k + 3;
end

if k <= n
    fprintf(out, '%.6f', data(k));
    k = k + 1;
    if k <= n
        fprintf(out, ' %.6f', data(k));
        k = k + 1;
    end
end

fprintf(out, '\n');
fclose(out);



