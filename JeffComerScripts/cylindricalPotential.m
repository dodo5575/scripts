% This Matlab script reads a .dx file and then plots
% averages the potential over the azimuthal angle.
% Also adds in an applied potential.
% Author: Jeff Comer <jcomer2@illinois.edu>
clear all;

%Parameters:
prefix = 'pore1.6_open_4V3-4';
%prefix = 'open_4V17-20';
ds = 0.5; % radial step in angstroms
dx = 1.0; % output grid interval in angstroms
gridX = 60.0; % size of X grid
thetaSamples = 90; % number of samples around the axis
potExternal = 4.0; % voltage drop along z in volts

% Read the file.
inName = strcat(prefix, '.dx');
outName = strcat(prefix, '_cyl.dx');
[data grid delta origin] = readdx(inName);
num = length(data);
if grid(1)*grid(2)*grid(3) ~= num
    error('cylindricalPotential:gridSize', 'Grid size does not match data size.')
end

% Get the bounds.
smax = 0.5*norm(grid(1)*delta(:,1) + grid(2)*delta(:,2));
zmin = origin(3);
ns = ceil(smax/ds);
nz = grid(3);
vs = ds*(0:ns-1)';
vz = zmin + delta(3,3)*(0:nz-1)';
[secZ secS] = meshgrid(vz,vs);

% Create the pot(s,z) map.
secV = zeros(ns, nz);
for it=1:thetaSamples
    ax = cos(2.0*it*pi/thetaSamples);
    ay = sin(2.0*it*pi/thetaSamples);
    for is=1:ns
        s = ds*(is-1);
        for iz=1:nz
            r = [s*ax; s*ay; zmin+(iz-1)*delta(3,3)];
            secV(is,iz) = secV(is,iz) +...
                interpolateLattice(data,grid,delta,origin,r);
        end
    end
end
secV = secV/thetaSamples;

% Parametrize an orthogonal lattice.
nx = ceil(gridX/dx);
dataOrtho = zeros(nx*nx*nz,1);
gridOrtho = [nx nx nz]';
deltaOrtho = [dx 0 0;0 dx 0; 0 0 delta(3,3)];
originOrtho = [-0.5*nx*dx -0.5*nx*dx origin(3)]';

% Compute the external potential in k T/e.
lz = delta(3,3)*nz;
external = -potExternal/0.0258522*vz/lz;

% Sample on a orthogonal lattice.
index = 1;
for ix=1:nx
    x = (ix-1)*dx + originOrtho(1);
    for iy=1:nx
        y = (iy-1)*dx + originOrtho(2);
        s = sqrt(x^2 + y^2);
        is = 1 + floor(s/ds);
        
        for iz=1:nz
            if is >= size(secV, 1)
                error('cylindricalPotential:outputGridSize',...
                    'Output grid size gridX is too large for input file')
            end
            s0 = secS(is, iz);
            s1 = secS(is+1, iz);
            f0 = secV(is, iz);
            f1 = secV(is+1, iz);
            a = (s-s0)/ds;
            
            dataOrtho(index) = f1*a + f0*(1.0-a) + external(iz);
            index = index + 1;
        end
    end
end

% Write the resulting dx file.
writedx(dataOrtho, gridOrtho, deltaOrtho, originOrtho, outName); 
dlmwrite(strcat(prefix,'_cyl.s'), secS, ' ');
dlmwrite(strcat(prefix,'_cyl.z'), secZ, ' ');
dlmwrite(strcat(prefix,'_cyl.pot'), secV, ' ');

% Plot the contours.
figure(1)
contourf(secS, secZ, secV, 14)
xlabel('s (A)')
ylabel('z (A)')
axis equal
colorbar



