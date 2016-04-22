% Perform a Monte Carlo simulation of a hairpin DNA.
% Author: Jeff Comer <jcomer2@illinois.edu>
clear all;

% Parameters:
%%%%%%%%%%%%%
% Simulation
nSteps = 100000;
movePeriod = 2;
outSuffix = 'str3';

%Display
progressPeriod = 100;
plotConform = 0;
dispBound = [-8 8 -8 8 -8 8];

% Dimensions
coilLength = 38.01; % nm
%helixLength = 3.529; % nm with loop
helixLength = 3.4;
ssDiam = 0.8; % nm
dsDiam = 2.2; % nm
% Zhang et al. Biophys. J. 81 (2001) 1133-1443.
kuhnLength = 1.6; % nm
strModulus = 123.5; % kB*T/nm^2

% Monte Carlo displacements
angleMax = pi;
distMax = 0.05*kuhnLength;

% Units and constants
kB = 1.3087e-23; % J/K
temp = 295; %K
qp = 1.6022e-19; % C
dielec = 78*8.85418782e-12; % F/m
%%%%%%%%%%%%%

% Initialize.
nNodes = ceil(coilLength/kuhnLength);
helixPos = [0 0 -helixLength; 0 0 0];
dist = 0.5*(ssDiam + dsDiam);
conformFile = sprintf('conform_ss%d_%s.dat', nNodes, outSuffix);
obscureFile = sprintf('obscure_ss%d_%s.dat', nNodes, outSuffix);
rmsdFile = sprintf('rmsd_ss%d_%s.dat', nNodes, outSuffix);
energyFile = sprintf('energy_ss%d_%s.dat', nNodes, outSuffix);
out = fopen(conformFile,'wt');

% Generate a random conformation.
reject = 1;
count = 0;
while reject
    r = randomConformation(nNodes, kuhnLength, 0.1);
    reject = excludeSingle(r,ssDiam);
    if ~reject
        reject = excludeDouble(helixPos,r,dist);
    end
    count = count + 1;
end
disp(sprintf('Attempts required to generate valid conformation: %d',count));

% Initialize statistics.
nObscure = zeros(nSteps,1);
bObscure = zeros(nSteps,1);
rmsd = zeros(nSteps,1);
energySamp = zeros(nSteps,1);
ri = r;

% Run the MC steps.
nRejected = 0;
energy = computeEnergyStretch(r,kuhnLength,strModulus);
for step=1:nSteps
    r0 = r;

    % Create a new conformation.
    if mod(step,movePeriod) ~= 0
        r = moveRotation(r,angleMax);
        r = moveStretch(r,distMax);
    else
        r = moveDoubletTriplet(r);
    end

    % Reject excluded moves.
    [reject collideA collideB] = excludeSingle(r,ssDiam);
    if ~reject
        [reject collideA collideB] = excludeDouble(helixPos,r,dist);
    end
    
    % Use the Metropolis criterion to reject moves.
    if ~reject
        energy0 = energy;
        energy = computeEnergyStretch(r,kuhnLength,strModulus);
        reject = computeMetropolis(energy-energy0);
        if reject
            energy = energy0;
        end
    end
    
    % Display progress.
    if mod(step,progressPeriod) == 0
        disp(sprintf('Progress: %f',step/nSteps));
        
        if plotConform
            figure(1)
            clf
            plot3(r(:,1), r(:,2), r(:,3), 'k.-',helixPos(:,1), helixPos(:,2), helixPos(:,3), 'b.-', ...
                collideA(:,1), collideA(:,2), collideA(:,3), 'ro-', collideB(:,1), collideB(:,2), collideB(:,3), 'ro-')
            axis(dispBound);
            axis vis3d
            xlabel('x (nm)')
            ylabel('y (nm)')
            zlabel('z (nm)')
            drawnow
        end
    end
    
    if reject
        nRejected = nRejected + 1;
        % Revert to old state.
        r = r0;
    end
    
    % Accumulate data.
    nObscure(step) = sum(r(:,3) < -helixLength);
    bObscure(step) = any(r(:,3) < -helixLength);
    rmsd(step) = sqrt(mean(sum((r-ri).^2,2)));
    energySamp(step) = energy;
    fprintf(out, '!%d\n', step-1);
    fprintf(out, '%.12g %.12g %.12g\n', r');
end
fclose(out);

% Display some results.
acceptanceRate = (nSteps-nRejected)/nSteps
obscureFraction = sum(bObscure/nSteps)
meanEnergy = mean(energySamp);
meanEnergyPerNode = meanEnergy/(nNodes-1)

figure(2)
t = (1:nSteps)';
plot(t,energySamp,'bo-')

dlmwrite(obscureFile, [t nObscure], ' ');
dlmwrite(rmsdFile, [t rmsd], ' ');
dlmwrite(energyFile, [t energySamp], ' ');



