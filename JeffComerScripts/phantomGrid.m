% phantomGrid.m
% Create a gridforce grid with a gaussian potential about each atom.
clear all;

coordFile = 'pore2.0_coords.txt';
outFile = 'pore2.0_nostick4.dx';
sigma = 4.0; % standard deviation of gaussian in angstroms
gridSpace = 1.5; % angstroms

% cutoff for pair interaction and buffer at top and bottom
cutoff = sigma*4.0;

% Load the atom coordinates.
rAtom = dlmread(coordFile, ' ');
nAtom = length(rAtom(:,1));

% Perform a cell decomposition on them.
[bag neigh origin cellGrid] = cellDecomposition(rAtom, cutoff);

% Compute the grid dimensions.
x0 = min(rAtom(:,1))-0.5*cutoff;
x1 = max(rAtom(:,1))+0.5*cutoff;
y0 = min(rAtom(:,2))-0.5*cutoff;
y1 = max(rAtom(:,2))+0.5*cutoff;
z0 = min(rAtom(:,3))-0.5*cutoff;
z1 = max(rAtom(:,3))+0.5*cutoff;

% Determine the number of grid points along the axes.
nx = floor((x1-x0)/gridSpace);
ny = floor((y1-y0)/gridSpace);
nz = floor((z1-z0)/gridSpace);
dx = (x1-x0)/nx;
dy = (y1-y0)/ny;
dz = (z1-z0)/nz;

n = nx*ny*nz;
disp(sprintf('Total grid nodes: %d', n));

% Determine the potential at each point on the grid.
alpha = 1.0/(2.0*sigma*sigma);
pot = zeros(n,1);
index = 1;
for ix=1:nx
    for iy=1:ny
        for iz=1:nz
            r = [x0+(ix-1)*dx y0+(iy-1)*dy z0+(iz-1)*dz];
            cellIndex = cellLookup(r, cutoff, origin, cellGrid);
            %disp(sprintf('%d %g %g %g', cellIndex, r(1), r(2), r(3)))
                        
            % Compute the potential using the points in all neighbors of
            % the cell.
            neighIndex = 1;
            curr = neigh(cellIndex, neighIndex);
            while curr > 0
                % Get the points in this cell.
                pos = bag{curr};
                
                if length(pos) > 0
                    % Compute the potential contribution
                    % at this point from the current cell.
                    dr = [pos(:,1)-r(1) pos(:,2)-r(2) pos(:,3)-r(3)];
                    drSq = dot(dr, dr, 2);
                    pot(index) = pot(index) + sum(exp(-alpha*drSq));
                end

                % Move to the next neighbor.
                neighIndex = neighIndex + 1;
                curr = neigh(cellIndex, neighIndex);
            end
            
            index = index + 1;
        end
    end
end

% Write a dx file based on this potential.
grid = [nx ny nz];
delta = [dx 0 0; 0 dy 0; 0 0 dz;];
origin = [x0 y0 z0];
writedx(pot, grid, delta, origin, outFile);



