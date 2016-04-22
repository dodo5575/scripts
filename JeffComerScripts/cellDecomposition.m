function [bag neigh origin grid] = cellDecomposition(pos, d)
% Organize a set of 3D points into cells.
% [BAG NEIGH ORIGIN GRID]= cellDecomposition(POS,D)
% Place a set of points (POS) in boxes of dimension D.
% Return the cell arrays of the points (bag) in each cell
% along with the domain of each cell as matrix
% and a cell array of the neighbors of each cell.
% Bags go x slow, y medium, z fast.
% Author: Jeff Comer <jcomer2@illinois.edu>

% Get the bounds.
xMin = min(pos(:,1));
xMax = max(pos(:,1));
yMin = min(pos(:,2));
yMax = max(pos(:,2));
zMin = min(pos(:,3));
zMax = max(pos(:,3));
origin = [xMin-d yMin-d zMin-d];

% Determine the number of cell arrays along each axis.
lx = ceil((xMax-xMin)/d)+2;
ly = ceil((yMax-yMin)/d)+2;
lz = ceil((zMax-zMin)/d)+2;
nBag = lx*ly*lz;
grid = [lx ly lz];

% Determine the neighbors.
neigh = zeros(nBag,28);
for j=1:nBag
    % Add each of the 27 neighbors if they exist.
    num = 1;
    for kx=-1:1
        for ky=-1:1
            for kz=-1:1
                index = j + kz + ky*lz + kx*lz*ly;
                if index >= 1 && index <= nBag
                    neigh(j,num) = index;
                    num = num + 1;
               end
            end
        end
    end
end

% Fill the bag with the points in each cell.
bag = cell(nBag,1);
n = length(pos(:,1));
for j=1:n
    index = cellLookup(pos(j,:), d, origin, grid);
    bag{index} = [bag{index};pos(j,:)];
end



