clear all;

% Parameters:
nSteps = 10000;
progressPeriod = 100;
outSuffix = 'run0';
coilLength = 38.01; % nm
%helixLength = 3.529; % nm with loop
helixLength = 3.4;
ssDiam = 0.8; % nm
dsDiam = 2.2; % nm
% Zhang et al. Biophys. J. 81 (2001) 1133-1443.
kuhnLength = 1.6; % nm

% Units and constants:
kB = 1.3087e-23; % J/K
temp = 295; %K
qp = 1.6022e-19; % C
dielec = 78*8.85418782e-12; % F/m

% Initialize.
nNodes = ceil(coilLength/kuhnLength);
helixPos = [0 0 -helixLength; 0 0 0];
dist = 0.5*(ssDiam + dsDiam);
conformFile = sprintf('conform_ss%1.1f_ds%1.1f_%s.dat', ssDiam, dsDiam, outSuffix);
obscureFile = sprintf('obscure_ss%1.1f_ds%1.1f_%s.dat', ssDiam, dsDiam, outSuffix);
rmsdFile = sprintf('rmsd_ss%1.1f_ds%1.1f_%s.dat', ssDiam, dsDiam, outSuffix);
out = fopen(conformFile,'wt');

% Generate a random conformation.
reject = 1;
count = 0;
while reject
    r = randomConformationRigid(nNodes, kuhnLength);
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
ri = r;

% Run the MC steps.
figure(1)
nRejected = 0;
for step=1:nSteps
    r0 = r;
    r = moveRotation(r,pi);
        
    [reject collideA collideB] = excludeSingle(r,ssDiam);
    if ~reject
        [reject collideA collideB] = excludeDouble(helixPos,r,dist);
    end
    
    % Display progress.
    if mod(step,progressPeriod) == 0
        disp(sprintf('Progress: %f',step/nSteps));
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
    fprintf(out, '!%d\n', step-1);
    fprintf(out, '%.12g %.12g %.12g\n', r');
end
acceptanceRate = (nSteps-nRejected)/nSteps
fclose(out);

plot3(r(:,1), r(:,2), r(:,3), 'k.-',helixPos(:,1), helixPos(:,2), helixPos(:,3), 'b.-', ...
            collideA(:,1), collideA(:,2), collideA(:,3), 'ro-', collideB(:,1), collideB(:,2), collideB(:,3), 'ro-')
axis vis3d
obscureFraction = sum(bObscure/nSteps)

figure(2)
t = (1:nSteps)';
plot(t,rmsd,'bo-')
dlmwrite(obscureFile, [t nObscure], ' ');
dlmwrite(rmsdFile, [t rmsd], ' ');




