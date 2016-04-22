% This Matlab script reads a .dx file and then interactively allows the
% user to plot different cross sections.
% Author: Jeff Comer <jcomer2@illinois.edu>
clear all;
% Input:
potFile = 'unzip0.dx';

[data grid delta origin] = readdx(potFile);
planeX = [1; 0; 0];
planeY = [0; 0; 1];
pos = [0; 0; 50];
planeNodes = 100;

% Convert the potential to V instead of k T/e.
data = 0.0258522*data;

plotting = 1;
currPlot = 1;
disp('Starting potential cross section plotter')
while plotting
disp('CTRL-C exits')
    
bad = 1;
while bad
    pos = input('Please give the position of the cross section plane (enclosed in []): ')
    if all(size(pos) == [1,3])
        pos = pos';
    end
    if all(size(pos) == [3,1])
        bad = 0;
    else
        disp('Bad input, please enter a vector')
    end
end

bad = 1;
while bad
    colat = input('Please give the polar angle of the plane (deg): ')
    if all(size(colat) == [1,1])
        bad = 0;
    else
        disp('Bad input, please enter a scalar')
    end
end

bad = 1;
while bad
    azi = input('Please give the azimuthal angle of the plane (deg): ')
    if all(size(azi) == [1,1])
        bad = 0;
    else 
        disp('Bad input, please enter a scalar')
    end
end


% Form the plane basis.
theta = pi/180.*colat;
phi = pi/180.*azi;
basis = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
basis = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]*basis;
planeX = basis(:,1);
planeY = basis(:,2);
planeZ = basis(:,3);
delta0 = delta;
delta0Inv = inv(delta0);

num = length(data);
if grid(1)*grid(2)*grid(3) ~= num
    error('plotpotential:gridSize', 'Grid size does not match data size.')
end

% Compute the extreme node positions in world space.
latticeR = zeros(3,7);
latticeR(:,1) = delta0*[grid(1); 0; 0] + origin;
latticeR(:,2) = delta0*[0; grid(2); 0] + origin;
latticeR(:,3) = delta0*[0; 0; grid(3)] + origin;
latticeR(:,4) = delta0*[0; grid(2); grid(3)] + origin;
latticeR(:,5) = delta0*[grid(1); 0; grid(3)] + origin;
latticeR(:,6) = delta0*[grid(1); grid(2); 0] + origin;
latticeR(:,7) = delta0*[grid(1); grid(2); grid(3)] + origin;
% Find the corresponding extrema in plane space.
planeMinPx = dot(planeX, origin-pos);
planeMinPy = dot(planeY, origin-pos);
planeMaxPx = dot(planeX, origin-pos);
planeMaxPy = dot(planeY, origin-pos);
for j=1:size(latticeR,2)
    px = dot(planeX, latticeR(:,j)-pos);
    py = dot(planeY, latticeR(:,j)-pos);
    if px < planeMinPx, planeMinPx = px; end
    if py < planeMinPy, planeMinPy = py; end
    if px > planeMaxPx, planeMaxPx = px; end
    if py > planeMaxPy, planeMaxPy = py; end
end


% Compute the cross section potentials.
secV = ones(planeNodes)*mean(data);
px = linspace(planeMinPx, planeMaxPx, planeNodes);
py = linspace(planeMinPy, planeMaxPy, planeNodes);
planeOrigin = planeMinPx*planeX + planeMinPy*planeY;
[secPy secPx] = meshgrid(py,px);
n = 1;
for ix=1:planeNodes
    for iy=1:planeNodes
        % Compute the plane node position in world space.
        r = pos + px(ix)*planeX + py(iy)*planeY;
        secV(ix,iy) = interpolateLattice(data,grid,delta,origin,r);
    end
end

% Plot the contours.
figure(currPlot)
contourf(secPx, secPy, secV, 14)
xlabel('plane x (A)')
ylabel('plane y (A)')
axis equal
colorbar
titleStr=sprintf('Potential at [%s %s %s] with polar angle %s and azimuth %s',...
    num2str(pos(1)), num2str(pos(2)), num2str(pos(3)), num2str(colat), num2str(azi));
title(titleStr)
currPlot = currPlot + 1;
end



