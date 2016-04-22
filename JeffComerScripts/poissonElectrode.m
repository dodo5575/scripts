% 
% jcomer2@uiuc.edu
clear all;

% Parameters:
h = 1.9; % in nm
plateX = 52;
plateY = 52;
outName = sprintf('electrode%d', plateX);
plateZ = 3*h;
plateDiam = 10;
poreZ = 20;
dnaZ = poreZ*1.5;
poreDiam = 4;
poreDiel = 5;
poreAngle = 10;
basepairLength = 0.34;
basepairCharge = -2; % in e
waterDiel = 80;

sysZ = 2.5*poreZ;
sysX = 2*plateX;
sysY = 2*plateY;
poreSlope = tan(10*pi/180);
totalCharge = dnaZ/basepairLength*basepairCharge;
chargeRad = 12*h/7;

% Units
eps0 = 8.854187817e-12;
elemCharge = 1.602176487e-19;

% Make the system.
outPrefix = sprintf('%s', outName);
[plan nx ny nz] = createSystem(sysX, sysY, sysZ, h);
x = plan(:,1);
y = plan(:,2);
z = plan(:,3);
ni = [nx ny nz];
nNodes = length(x);
disp(sprintf('Created a system of %d nodes.', nNodes));

% Make the dielectric map.
disp('Making the dielectric map.');
eps = zeros(size(x));
for k=1:nNodes
  s = 0.5*poreDiam + poreSlope*abs(z(k));
  if abs(z(k)) < 0.5*poreZ && x(k)^2 + y(k)^2 > s^2 
    eps(k) = poreDiel;
  else
    eps(k) = waterDiel;
  end
end

% Get the gradient of the dielectric map.
disp('Computing the gradient of the dielectric map.');
eps = arrayBlur3d(eps, nx, ny, nz);
eps = arrayBlur3d(eps, nx, ny, nz);
[epsX epsY epsZ] = arrayGrad3d(eps, nx, ny, nz, h);

% Make a conductance map.
disp('Making the conductance map...');
condNode = zeros(size(x));
condCount = 0;
cond = zeros(size(x));
for k=1:nNodes
  if abs(x(k)) < 0.5*plateX && abs(y(k)) < 0.5*plateY && abs(z(k)) < 0.5*plateZ
    if x(k)^2 + y(k)^2 < 0.25*plateDiam^2, continue, end
    
    cond(k) = 1;
    condCount = condCount + 1;
    condNode(condCount) = k;
  end
end
condNode = condNode(1:condCount);
%cond = zeros(size(x));
disp(sprintf('Found %d conductor nodes.', condCount));


% Make the charge map.
disp('Making the charge density map...');
rho = zeros(size(x));
nCharges = 0;
for j=1:nNodes
  if x(j)^2 + y(j)^2 < chargeRad^2 && abs(z(j)) < 0.5*dnaZ
    rho(j) = 1;
    nCharges = nCharges + 1;
  end
end
rho = rho*totalCharge/(nCharges*h^3); % in e nm^-3
qStar = elemCharge/eps0*1e9/waterDiel;

% Form the relation matrix m and the vector b in mv = b;
m = sparse(nNodes,nNodes);
b = zeros(nNodes,1);

disp(sprintf('\nFilling the matrix...'))
disp('Setting the surface boundary conditions...')
[m b] = zeroBoundary(m, b, nx, ny, nz);

disp('Setting the internal nodes...')
surfCount = 0;
surfNode = zeros(size(x));
insideCount = 0;
insideNode = zeros(size(x));
% Fill the nodes not on the edges.
l = 0;
for iz=0:nz-1
  for iy=1:ny-2
    for ix=1:nx-2
      l = l + 1;
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
      m(j,k(2)) = 1 - 0.5*h*epsX(j)/eps(j);
      m(j,k(3)) = 1 + 0.5*h*epsX(j)/eps(j); 
      m(j,k(4)) = 1 - 0.5*h*epsY(j)/eps(j); 
      m(j,k(5)) = 1 + 0.5*h*epsY(j)/eps(j); 
      m(j,k(6)) = 1 - 0.5*h*epsZ(j)/eps(j); 
      m(j,k(7)) = 1 + 0.5*h*epsZ(j)/eps(j); 
      
      b(j) = -h^2*qStar*rho(j);
    end
  end
  disp(sprintf('%d percent complete', floor(100*l/nNodes)));
end
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

% Solve the system.
disp('Solving Poisson`s equation');
v = m\b;

% Find the average value at the conducting node.
condV = mean(v(condNode))
condVStd = std(v(condNode))

% Compute the gradient.
[ex ey ez] = arrayGrad3dNonperiodic(v, nx, ny, nz, h);
ex = -ex;
ey = -ey;
ez = -ez;

% Extract a cross section.
secY = zeros(ny,nz);
secZ = zeros(ny,nz);
secV = zeros(ny,nz);
secEy = zeros(ny,nz);
secEz = zeros(ny,nz);
ixc = floor(nx/2);
for iz=0:nz-1
  for iy=0:ny-1
    j = 1 + ixc + nx*iy + nx*ny*iz;
    secY(iy+1, iz+1) = y(j);
    secZ(iy+1, iz+1) = z(j);
    secV(iy+1, iz+1) = v(j);
    secEy(iy+1, iz+1) = ey(j);
    secEz(iy+1, iz+1) = ez(j);
  end
end

% Extract a profile.
izc = floor(nz/2);
iyc = floor(ny/2);
profileZ = zeros(nx,1);
profileV = zeros(nx,1);
profileEy = zeros(nx,1);
profileEz = zeros(nx,1);
for ix=0:nx-1
  j = 1 + ix + nx*iyc + nx*ny*izc;
  profileZ(ix+1) = x(j);
  profileV(ix+1) = v(j);
  profileEy(ix+1) = ey(j);
  profileEz(ix+1) = ex(j);
end

% Write the results.
disp('Writing the results...')
dlmwrite(strcat(outPrefix,'.y.sec'), secY, ' ');
dlmwrite(strcat(outPrefix,'.z.sec'), secZ, ' ');
dlmwrite(strcat(outPrefix,'.v.sec'), secV, ' ');
dlmwrite(strcat(outPrefix,'.ey.sec'), secEy, ' ');
dlmwrite(strcat(outPrefix,'.ez.sec'), secEz, ' ');
dlmwrite(strcat(outPrefix,'.v.dat'), [profileZ profileV], ' ');
dlmwrite(strcat(outPrefix,'.ey.dat'), [profileZ profileEy], ' ');
dlmwrite(strcat(outPrefix,'.ez.dat'), [profileZ profileEz], ' ');

% Make dx file.
% Rearrage the data to z fast, y medium, x slow.
disp('Writing the dx file.')
dxData = zeros(size(v));
n = 0;
for ix=0:nx-1
  for iy=0:ny-1
    for iz=0:nz-1
      j = 1 + ix + nx*iy + nx*ny*iz;
      n = n + 1;
      dxData(n) = v(j);
    end
  end
end
writedx(dxData, [nx ny nz], h*eye(3), [x(1) y(1) z(1)], strcat(outPrefix,'.dx'));


% Plot things.
figure(2)
clf

%contourf(secY, secZ, secV0, linspace(0,0.1,10))
contourf(secY, secZ, secV, 14)
colormap(jet)
xlabel('y (nm)')
ylabel('z (nm)')
axis square
colorbar

figure(3)
clf
plot(profileZ, profileV, 'b--')
xlabel('z (nm)')
ylabel('u (V)')
%axis([-inf inf -inf inf])

figure(1)
clf
hold on
plot3(x(insideNode), y(insideNode), z(insideNode), 'ko')
plot3(x(surfNode), y(surfNode), z(surfNode), 'bo')
hold off
axis vis3d