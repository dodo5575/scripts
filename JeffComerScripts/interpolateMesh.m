function n = interpolateMesh(gridSpace, xFile, yFile, zFile, potFile, outFile)
% Interpolate a field from a grid onto a lattice and write a dx file.
% The input grid vectors are orthogonal but not necessarily equally-spaced.
% The potential data file must be structured x fast, y medium, z slow.
% GRIDSPACE is the side length of the sampling cubes.

% Read the data.
x = dlmread(xFile);
y = dlmread(yFile);
z = dlmread(zFile);
pot = dlmread(potFile);
mesh = [length(x) length(y) length(z)];

% Calculate the resulting lattice dimensions.
origin = [x(1)+0.5*gridSpace y(1)+0.5*gridSpace z(1)+0.5*gridSpace];
nx = floor((x(end) - x(1))/gridSpace)-1;
ny = floor((y(end) - y(1))/gridSpace)-1;
nz = floor((z(end) - z(1))/gridSpace)-1;
n = nx*ny*nz;

% Sample the mesh.
dxPot = zeros(n,1);
j = 1;
for ix=1:nx
    for iy=1:ny
        for iz=1:nz
            % Compute the position of the lattice node.
            xj = ix*gridSpace + origin(1);
            yj = iy*gridSpace + origin(2);
            zj = iz*gridSpace + origin(3);
            
            % Determine the mesh vertices that neighbor this node.
            % This is faster with a binary search instead of a linear
            % search.
            % Make the p values zero-based.
            px1 = binsearch(xj,x);
            px0 = px1-1;
            py1 = binsearch(yj,y);
            py0 = py1-1;
            pz1 = binsearch(zj,z);
            pz0 = pz1-1;
            
            % Compute the interpolation coefficients.
            a = (xj-x(px0+1))/(x(px1+1)-x(px0+1));
            b = (yj-y(py0+1))/(y(py1+1)-y(py0+1));
            c = (zj-z(pz0+1))/(z(pz1+1)-z(pz0+1));
                   
            % Extract the values of the potential on the neighboring vertices.
            f000 = pot(1 + px0 + mesh(3)*py0 + mesh(2)*mesh(3)*pz0);
            f001 = pot(1 + px1 + mesh(3)*py0 + mesh(2)*mesh(3)*pz0);
            f010 = pot(1 + px0 + mesh(3)*py1 + mesh(2)*mesh(3)*pz0);
            f011 = pot(1 + px1 + mesh(3)*py1 + mesh(2)*mesh(3)*pz0);
            f100 = pot(1 + px0 + mesh(3)*py0 + mesh(2)*mesh(3)*pz1);
            f101 = pot(1 + px1 + mesh(3)*py0 + mesh(2)*mesh(3)*pz1);
            f110 = pot(1 + px0 + mesh(3)*py1 + mesh(2)*mesh(3)*pz1);
            f111 = pot(1 + px1 + mesh(3)*py1 + mesh(2)*mesh(3)*pz1);
            
            % Mix along lattice vector x.
            g00 = f100*a + (1.0-a)*f000;
            g10 = f110*a + (1.0-a)*f010;
            g01 = f101*a + (1.0-a)*f001;
            g11 = f111*a + (1.0-a)*f011;

            % Mix along lattice vector y.
            h0 = g10*b + (1.0-b)*g00;
            h1 = g11*b + (1.0-b)*g01;

            % Mix along lattice vector z.
            val = h1*c + (1.0-c)*h0;
            dxPot(j) = val;
            
            j = j + 1;
        end
    end
end

writedx(dxPot, [nx ny nz], gridSpace*eye(3), origin, outFile);




