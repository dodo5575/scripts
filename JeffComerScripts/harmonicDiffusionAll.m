clear all;

indexList = 0:41;
sys = 'at_1M_pot';
indexFile = sprintf('diff_window_%s.txt', sys);
dt = 16e-6; % in ns
subSteps = 600; % steps in a subtrajectory
cutIndex = 10;
inDir = 'traj';
diffusionCut = 1000;

fprintf('\n****%s\n', datestr(now));
inPrefix = sprintf('%s/bpDiff_%s_pos', inDir, sys);
outPrefix = sprintf('try/hummer_%s_pos', sys);

diffuse = zeros(length(indexList),3);
for ind=1:length(indexList)
    index = indexList(ind);
    
    inFileGlob = sprintf('%s%d-*.dat', inPrefix, index);
    outPre = sprintf('%s%d', outPrefix, index);
    
    % Read the reference position.
    refPos = readRefPosition(indexFile, index);
    
    % Get the list of files.
    inFileList = dir(inFileGlob);
    %inFileList = inFileList(1:2);
    
    % Load the data.
    posList = [];
    deltaCoorX = [];
    deltaCoorY = [];
    deltaCoorZ = [];
    fprintf('\n');
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
        
        if isempty(data)
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
    
    if isempty(posList)
        continue
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
    diffusionX = diffusionX(diffusionCut:end);
    diffTimeX = tauTimeX(diffusionCut:end);
   
    diffusionY = varY./tauValY;
    diffusionY = diffusionY(diffusionCut:end);
    diffTimeY = tauTimeY(diffusionCut:end);
    
    diffusionZ = varZ./tauValZ;
    diffusionZ = diffusionZ(diffusionCut:end);
    diffTimeZ = tauTimeZ(diffusionCut:end);
    
    % Write the results.
    fprintf('INDEX: %d\n', index);
    fprintf('SUBTIME: %d\n', subSteps*dt);
    fprintf('SUBSTEPS: %d\n', subSteps);
    fprintf('SAMPLES: %d\n', length(deltaCoorX));
    fprintf('TOTALTIME: %d\n', length(deltaCoorX)*subSteps*dt);
    fprintf('POINTS: %d\n', length(posList));
    fprintf('VAR: %g %g %g\n', varX, varY, varZ);
    fprintf('STD: %g %g %g\n', sqrt(varX), sqrt(varY), sqrt(varZ));
    fprintf('MEAN: %g %g %g\n', posX, posY, posZ);
    fprintf('REF: %g %g %g\n', refPos(1), refPos(2), refPos(3));
    fprintf('DIFFUSE: %g %g %g\n', diffusionX(end), diffusionY(end), diffusionZ(end));
    
    out = fopen(strcat(outPre,'.log'), 'w');
    fprintf(out, 'INDEX: %d\n', index);
    fprintf(out, 'SUBTIME: %d\n', subSteps*dt);
    fprintf(out, 'SUBSTEPS: %d\n', subSteps);
    fprintf(out, 'SAMPLES: %d\n', length(deltaCoorX));
    fprintf(out, 'TOTALTIME: %d\n', length(deltaCoorX)*dt);
    fprintf(out, 'POINTS: %d\n', length(posList));
    fprintf(out, 'VAR: %g %g %g\n', varX, varY, varZ);
    fprintf(out, 'STD: %g %g %g\n', sqrt(varX), sqrt(varY), sqrt(varZ));
    fprintf(out, 'MEAN: %g %g %g\n', posX, posY, posZ);
    fprintf(out, 'REF: %g %g %g\n', refPos(1), refPos(2), refPos(3));
    fprintf(out, 'DIFFUSE: %g %g %g\n', diffusionX(end), diffusionY(end), diffusionZ(end));
    fclose(out);
    
    diffuse(ind,1) = diffusionX(end);
    diffuse(ind,2) = diffusionY(end);
    diffuse(ind,3) = diffusionZ(end);
        
    % Write output.
    dlmwrite(strcat(outPre,'_corrNormX.dat'), [tim corrNormX], ' ');
    dlmwrite(strcat(outPre,'_corrNormY.dat'), [tim corrNormY], ' ');
    dlmwrite(strcat(outPre,'_corrNormZ.dat'), [tim corrNormZ], ' ');
    
    dlmwrite(strcat(outPre,'_tauValX.dat'), [tauTimeX tauValX], ' ');
    dlmwrite(strcat(outPre,'_tauValY.dat'), [tauTimeY tauValY], ' ');
    dlmwrite(strcat(outPre,'_tauValZ.dat'), [tauTimeZ tauValZ], ' ');
    
    dlmwrite(strcat(outPre,'_diffusionX.dat'), [diffTimeX diffusionX], ' ');
    dlmwrite(strcat(outPre,'_diffusionY.dat'), [diffTimeY diffusionY], ' ');
    dlmwrite(strcat(outPre,'_diffusionZ.dat'), [diffTimeZ diffusionZ], ' ');
end

dlmwrite(strcat(outPrefix,'_diffuse.dat'), [indexList' diffuse], ' ');
