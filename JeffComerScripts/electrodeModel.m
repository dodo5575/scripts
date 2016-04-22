% electrodeModel.m
clear all;

probeRadius = 21.882; % in A
conc = 0.1; % in mol/kg
dl = 0.5; % in A

srcLambda0 = 2/3.4; % DNA charge density in e/A
srcDiameter = 20; % DNA diameter in A
srcLength = 200; % in A

% Constants
elemCharge = 1.602176487e-19; % in C
eps0 = 8.854187817e-12; % in SI
debye = 0.3163*sqrt(conc); % in A^-1
eps = 80.0;
elecConst = elemCharge/(4*pi*eps0*eps);

boltzmann = 1.308e-23; % J/K
avogadro = 6.02e23; 
density = conc*avogadro*1000; % m^-3
debyeKappa = elemCharge/sqrt(boltzmann*300*eps*eps0/density);
kappaConst = 1e-10*elemCharge/sqrt(boltzmann*300*eps*eps0/avogadro/1000);

molSize = 2:0.1:5;
outLen = zeros(length(molSize),1);
outPot = zeros(length(molSize),1);
for j=1:length(molSize)
srcLambda = -2/molSize(j); % DNA in e/A
    
% Generate a cross section of the src cylinder.
nx = ceil(srcDiameter/dl);
ny = nx;
srcRadSq = (0.5*srcDiameter)^2;
originX = -0.5*srcDiameter;
originY = originX;
srcX = zeros(nx*ny,1);
srcY = zeros(nx*ny,1);
srcN = 0;
for ix=1:nx
    for iy=1:ny
        x = originX + ix*dl;
        y = originY + iy*dl;
        
        if x*x + y*y < srcRadSq
            srcN = srcN + 1;
            srcX(srcN) = x;
            srcY(srcN) = y;
        end
    end
end

% Reform the src arrays.
srcX = srcX(1:srcN);
srcY = srcY(1:srcN);

originZ = -0.5*srcLength;
nz = ceil(srcLength/dl);

%figure(1)
%plot(srcX, srcY, '.')

pot = 0;
probeX = probeRadius;
for iz=1:nz
    for ia=1:srcN
        d = sqrt((probeX-srcX(ia))^2 + srcY(ia)^2 + (originZ + iz*dl)^2);
        pot = pot + exp(-debye*d)/d;
    end
end
pot = elecConst*pot*dl*srcLambda/srcN*1e10*1000;
%sprintf('size: %f, %f mV', molSize, pot)
outLen(j) = molSize(j);
outPot(j) = pot;
end

%plot(1./outLen, outPot)
invLen = 1./(outLen*0.1); % in nm^-1
plot(invLen, outPot)
slope = (outPot(end)-outPot(1))/(invLen(end)-invLen(1))



