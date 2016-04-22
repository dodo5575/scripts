clear all;

index = 40;
indexFile = 'diff_window.txt';
inDir = 'traj';
sys = 'at_chl';
inFileGlob = sprintf('%s/bpDiff_%s_pos%d-*.dat', inDir, sys, index);
outPrefix = sprintf('hummer/hummer_%s_pos%d.dat', sys, index);
dt = 16e-6; % in ns
subSteps = 600; % steps in a subtrajectory
cutIndex = 10;

fprintf('\n****%s\n', datestr(now));
deltaCoorX = [];
deltaCoorY = [];
deltaCoorZ = [];
subNum = 0;

% Read the reference position.
refPos = readRefPosition(indexFile, index);

% Get the list of files.
inFileList = dir(inFileGlob);
%inFileList = inFileList(1:2);

% Load the data.
posList = [];
for k=1:length(inFileList)
    inFile = sprintf('%s/%s', inDir, inFileList(k).name);
    fprintf('Loading %s...\n', inFile);
    data = dlmread(inFile, ' ');
    
    % Trim the data.
    data = data(cutIndex:end,:);
    n = length(data(:,1));
    fprintf('Loaded %d samples.\n', n);
    posList = [posList; data];
    
    subCount = floor(n/subSteps);
    subStepCount = subCount*subSteps;
    
    if length(data) < 1
        continue
    end
    
    % Chop into subtrajectories.
    subX = reshape(data(1:subStepCount,1), subSteps, subCount)';
    subY = reshape(data(1:subStepCount,2), subSteps, subCount)';
    subZ = reshape(data(1:subStepCount,3), subSteps, subCount)';
    
    deltaCoorX = [deltaCoorX; subX];
    deltaCoorY = [deltaCoorY; subY];
    deltaCoorZ = [deltaCoorZ; subZ];
end

% Compute the variance and means.
posX = mean(posList(:,1));
posY = mean(posList(:,2));
posZ = mean(posList(:,3));
varX = var(posList(:,1)-posX);
varY = var(posList(:,2)-posY);
varZ = var(posList(:,3)-posZ);

% Shift by the mean.
deltaCoorX = deltaCoorX - posX;
deltaCoorY = deltaCoorY - posY;
deltaCoorZ = deltaCoorZ - posZ;

%for k=1:length(deltaCoorX(:,1))
%    deltaCoorX() = deltaCoorX - mean(
%end

fprintf('INDEX: %d\n', index);
fprintf('SUBTIME: %d\n', subSteps*dt);
fprintf('SUBSTEPS: %d\n', subSteps);
fprintf('SAMPLES: %d\n', length(deltaCoorX));
fprintf('POINTS: %d\n', length(posList));
fprintf('VAR: %g %g %g\n', varX, varY, varZ);
fprintf('STD: %g %g %g\n', sqrt(varX), sqrt(varY), sqrt(varZ));
fprintf('MEAN: %g %g %g\n', posX, posY, posZ);
fprintf('REF: %g %g %g\n', refPos(1), refPos(2), refPos(3));

% Compute the mean correlation function.
tim = (0:(subSteps-1))'*dt;
corrX = zeros(subSteps,1);
corrY = zeros(subSteps,1);
corrZ = zeros(subSteps,1);
% Loop over each timestep.
for s=1:subSteps
    corrX(s) = mean(deltaCoorX(:,1).*deltaCoorX(:,s));
    corrY(s) = mean(deltaCoorY(:,1).*deltaCoorY(:,s));
    corrZ(s) = mean(deltaCoorZ(:,1).*deltaCoorZ(:,s));
end

corrNormX = corrX/varX;
corrNormY = corrY/varY;
corrNormZ = corrZ/varZ;

[tauTimeX tauValX] = integrateCorrelation(tim, corrNormX, 20);
[tauTimeY tauValY] = integrateCorrelation(tim, corrNormY, 20);
[tauTimeZ tauValZ] = integrateCorrelation(tim, corrNormZ, 20);

diffusionX = varX./tauValX;
diffusionY = varY./tauValY;
diffusionZ = varZ./tauValZ;

fprintf('DIFFUSE: %g %g %g\n', diffusionX(end), diffusionY(end), diffusionZ(end));

% Write output.
dlmwrite(strcat(outPrefix,'_corrNormX.dat'), [tim corrNormX], ' ');
dlmwrite(strcat(outPrefix,'_corrNormY.dat'), [tim corrNormY], ' ');
dlmwrite(strcat(outPrefix,'_corrNormZ.dat'), [tim corrNormZ], ' ');

dlmwrite(strcat(outPrefix,'_tauValX.dat'), [tauTimeX tauValX], ' ');
dlmwrite(strcat(outPrefix,'_tauValY.dat'), [tauTimeY tauValY], ' ');
dlmwrite(strcat(outPrefix,'_tauValZ.dat'), [tauTimeZ tauValZ], ' ');

dlmwrite(strcat(outPrefix,'_diffusionX.dat'), [tauTimeX diffusionX], ' ');
dlmwrite(strcat(outPrefix,'_diffusionY.dat'), [tauTimeY diffusionY], ' ');
dlmwrite(strcat(outPrefix,'_diffusionZ.dat'), [tauTimeZ diffusionZ], ' ');

% Plot.
figure(1)
plot(tim, corrNormX, 'r', tim, corrNormY, 'g', tim, corrNormZ, 'b');
legend('x', 'y', 'z');
xlabel('time (ns)');
ylabel('Q_{corr}/var(Q)')

figure(2)
plot(tauTimeX, tauValX, 'r', tauTimeY, tauValY, 'g', tauTimeZ, tauValZ, 'b');
legend('x', 'y', 'z');
xlabel('time (ns)');
ylabel('correlation time tau (ns)')

figure(3)
st = 1000;
plot(tim, corrNormX, 'r', tim, corrNormY, 'g', tim, corrNormZ, 'b');
plot(tauTimeX(st:end), varX./tauValX(st:end), 'r', tauTimeY(st:end), varY./tauValY(st:end), 'g', tauTimeZ(st:end), varZ./tauValZ(st:end), 'b');
legend('x', 'y', 'z');
xlabel('time (ns)');
ylabel('diffusivity (A^2/ns)')
