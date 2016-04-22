% This Matlab script reads a .dx file.
% Author: Jeff Comer <jcomer2@illinois.edu>
clear all;

%Parameters:
prefix = 'open_6V6-9';

% Read the file.
inName = strcat(prefix, '.dx');
outName = strcat(prefix, '_ortho.dx');
[data grid delta origin] = readdx(inName);
num = length(data);
if grid(1)*grid(2)*grid(3) ~= num
    error('orthogonalPotential:gridSize', 'Grid size does not match data size.')
end

% Get the bounds.
l = grid(1)*delta(:,1) + grid(2)*delta(:,2) + grid(3)*delta(:,3);
gridOrtho = grid;
originOrtho = origin;
deltaOrtho = [l(1) 0 0; 0 l(2) 0; 0 0 l(3)];
dataOrtho = zeros(grid(1)*grid(2)*grid(3),1);

% Interpolate the potential on an orthogonal lattice.
disp('Interpolating onto orthogonal lattice...');
index = 1;
for ix=0:grid(1)
    for iy=0:grid(2)
        for iz=0:grid(3)
            r = originOrtho + (ix-1)*l(1) + (iy-1)*l(2) + (iz-1)*l(3);
            dataOrtho(index) = interpolateLattice(data,grid,delta,origin,r);
            index = index + 1;
        end
    end
end
disp('Finished with interpolation.')

% Write the resulting dx file.
writedx(dataOrtho, gridOrtho, deltaOrtho, originOrtho, outName); 




