% 
% jcomer2@uiuc.edu
clear all;

% Parameters:
outName = 'cond_sphere';
h = 0.24; % in nm
sphereZ = -0.1;
sphereRad = 1; % in nm
side = 5.5;
chargeZ = 1.7;
chargeVal = 1; % in e
chargeRad = 9*h/7;

% Compute things.
dist = abs(chargeZ - sphereZ);
outPrefix = sprintf('%s_%d', outName, sphereRad);
[plan nx ny nz] = createSystem(side, side, 2*side, h);
x = plan(:,1);
y = plan(:,2);
z = plan(:,3);
ni = [nx ny nz];
nNodes = length(x);
disp(sprintf('Created a system of %d nodes.', nNodes));
disp(sprintf('Distance between charge and sphere: %f.', dist));

% Make a conductance map.
condNode = zeros(size(x));
condCount = 0;
cond = zeros(size(x));
for j=1:nNodes
  if x(j)^2 + y(j)^2 + (z(j)-sphereZ)^2 < (sphereRad+h)^2
    cond(j) = 1;
    condCount = condCount + 1;
    condNode(condCount) = j;
  end
end
condNode = condNode(1:condCount);
disp(sprintf('Found %d conductor nodes.\n', condCount));

% Units
eps0 = 8.854187817e-12;
elemCharge = 1.602176487e-19;

% Make the charge map.
rho = zeros(size(x));
nCharges = 0;
for j=1:nNodes
  if x(j)^2 + y(j)^2 + (z(j)-chargeZ)^2 < chargeRad^2
    rho(j) = 1;
    nCharges = nCharges + 1;
  end
end
rho = rho*chargeVal/(nCharges*h^3); % in e nm^-3
qStar = elemCharge/eps0*1e9;

% Form the relation matrix m and the vector b in mv = b;
m = sparse(nNodes,nNodes);
b = zeros(nNodes,1);

disp('Filling the matrix...')
disp('Setting the surface boundary conditions...')
[m b] = zeroBoundary(m, b, nx, ny, nz);

figure(1)
clf
axis vis3d
hold on

disp('Setting the internal nodes...')
surfCount = 0;
surfNode = zeros(size(x));
insideCount = 0;
insideNode = zeros(size(x));
% Fill the nodes not on the edges
for iz=1:nz-2
  for iy=1:ny-2
    for ix=1:nx-2
      % Add the contributions to the matrix.
      k = getNeighbors([ix iy iz], ni);
      j = k(1);
      
      inside = sum(cond(k));    
      if cond(j)
	if inside < 7
	  surfCount = surfCount + 1;
	  surfNode(surfCount) = j;
	else
	  insideCount = insideCount + 1;
	  insideNode(insideCount) = j;
	end
      end
      
      % Insulator node or internal conducting node.
      m(j,j) = -6;
      m(j,k(2)) = 1;
      m(j,k(3)) = 1; 
      m(j,k(4)) = 1; 
      m(j,k(5)) = 1; 
      m(j,k(6)) = 1; 
      m(j,k(7)) = 1; 
      
      b(j) = -h^2*qStar*rho(j);
    end
  end
end
hold off
surfNode = surfNode(1:surfCount);
insideNode = insideNode(1:insideCount);
disp(sprintf('Found %d conductor surface nodes.', surfCount));
disp(sprintf('Found %d conductor internal nodes.', insideCount));

disp('Setting the conductor internal nodes...')
% Make the internal nodes the average of the surface nodes.
l = insideNode(1);
b(l) = 0;
m(l,:) = 0;
m(l,surfNode) = 1;
m(l,l) = -(surfCount-1);
l0 = l;
for j=2:insideCount
  l = insideNode(j);
  b(l) = 0;
  m(l,:) = 0;
  m(l,l0) = 1;
  m(l,l) = -1;
  
  l0 = l;
end

% The whole matrix should be filled. Check to be sure.
if ~all(any(m))
  disp(sprintf('The matrix is missing entries!'));
end

% Set the analytic solution.
disp('Calculating the analytic solution...')
r = sqrt(x.^2+y.^2+(z-sphereZ).^2);
cost = (z-sphereZ)./r;
d = dist;
R = sphereRad;
%v0 = qStar/(4*pi)./sqrt(x.^2 + y.^2 + (z-chargeZ).^2);
v0 = qStar/(4*pi)*(1./sqrt(r.^2 + d^2 - 2*d*r.*cost) - 1./ ...
			sqrt(R^2 + d^2/R^2*r.^2 - 2*d*r.*cost) + ...
		   R./(d*r));
vR = qStar/(4*pi)/d;
v0 = (r <= R)*vR + (r > R).*v0;

% Solve the system.
disp('Solving Poisson`s equation');
v = m\b;

% Find the average value at the conducting node.
condV = mean(v(condNode))
condVStd = std(v(condNode))

%v(surfNode) = max(v(surfNode));

% Compute the gradient.
[ex ey ez] = arrayGrad3dNonperiodic(v, nx, ny, nz, h);
ex = -ex;
ey = -ey;
ez = -ez;

% Extract a cross section.
secY = zeros(ny,nz);
secZ = zeros(ny,nz);
secV = zeros(ny,nz);
secV0 = zeros(ny,nz);
secEy = zeros(ny,nz);
secEz = zeros(ny,nz);
ixc = floor(nx/2);
for iz=0:nz-1
  for iy=0:ny-1
    j = 1 + ixc + nx*iy + nx*ny*iz;
    secY(iy+1, iz+1) = y(j);
    secZ(iy+1, iz+1) = z(j);
    secV(iy+1, iz+1) = v(j);
    secV0(iy+1, iz+1) = v0(j);
    secEy(iy+1, iz+1) = ey(j);
    secEz(iy+1, iz+1) = ez(j);
  end
end

% Extract a profile.
ixc = floor(nx/2);
iyc = floor(ny/2);
profileZ = zeros(nz,1);
profileV = zeros(nz,1);
profileV0 = zeros(nz,1);
profileEy = zeros(nz,1);
profileEz = zeros(nz,1);
for iz=0:nz-1
  j = 1 + ixc + nx*iyc + nx*ny*iz;
  profileZ(iz+1) = z(j);
  profileV(iz+1) = v(j);
  profileV0(iz+1) = v0(j);
  profileEy(iz+1) = ey(j);
  profileEz(iz+1) = ez(j);
end

% Write the results.
disp('Writing the results...')
dlmwrite(strcat(outPrefix,'.y.sec'), secY, ' ');
dlmwrite(strcat(outPrefix,'.z.sec'), secZ, ' ');
dlmwrite(strcat(outPrefix,'.v.sec'), secV, ' ');
dlmwrite(strcat(outPrefix,'.v0.sec'), secV0, ' ');
dlmwrite(strcat(outPrefix,'.ey.sec'), secEy, ' ');
dlmwrite(strcat(outPrefix,'.ez.sec'), secEz, ' ');
dlmwrite(strcat(outPrefix,'.v.dat'), [profileZ profileV], ' ');
dlmwrite(strcat(outPrefix,'.v0.dat'), [profileZ profileV0], ' ');
dlmwrite(strcat(outPrefix,'.ey.dat'), [profileZ profileEy], ' ');
dlmwrite(strcat(outPrefix,'.ez.dat'), [profileZ profileEz], ' ');

% Make dx file.
% Rearrage the data to z fast, y medium, x slow.
disp('Writing the dx file.')
dxData = zeros(size(v));
dxData0 = zeros(size(v0));
n = 0;
for ix=0:nx-1
  for iy=0:ny-1
    for iz=0:nz-1
      j = 1 + ix + nx*iy + nx*ny*iz;
      n = n + 1;
      dxData(n) = v(j);
      dxData0(n) = v0(j);
    end
  end
end
writedx(dxData, [nx ny nz], h*eye(3), [x(1) y(1) z(1)], strcat(outPrefix,'.dx'));
writedx(dxData0, [nx ny nz], h*eye(3), [x(1) y(1) z(1)], ...
	strcat(outPrefix,'.0.dx'));

% Plot things.
figure(2)
clf

%contourf(secY, secZ, secV0, linspace(0,0.1,10))
contourf(secY, secZ, secV, 14)
colormap(jet)
xlabel('y (nm)')
ylabel('z (nm)')
axis equal
colorbar

figure(3)
clf
plot(profileZ, profileV, 'b--', profileZ, profileV0, 'k')
xlabel('z (nm)')
ylabel('u (V)')
axis([-inf inf 0 max(profileV)])

%figure(3)
%cen = floor(nx/2) + nx*floor(ny/2);
%ind = cen:(nx*ny):nNodes-1;
%plot(z(ind), v(ind), 'bo-')
%xlabel('z (nm)')
%ylabel('v (nm)')
