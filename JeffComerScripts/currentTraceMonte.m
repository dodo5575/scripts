clear all;

inDir = 'output';
inPrefix = 'curr_coil_alek_6V1-5';
voltage = 4;
outDir = 'trace';
nBins = 10;
step = 2;
expTempInitial = -1;
expTempFinal = -8;
nAttempts = 20000;

current0 = [];
%current0 = [1.0 1.45 0.5; 1.5 2.18 0.5; 2.0 2.91 0.5; 2.5 3.63 0.5; 3.0 4.36 0.5; 4.0 5.81 0.5; 5.0 7.26 0.5; 6.0 8.72 0.5; 8.0 11.6 0.5]; % linear
%current0 = [1.0 2.7 0.3; 1.5 3.9 0.1; 2.0 4.8 0.3; 2.5 5.0 1.0; ...
    %3.0 5.18 0.3; 4.0 5.81 0.08; 5.0 6.7 1.0; 6.0 7.5 0.1; 8.0 10.0 1.0];
currFactor = 1.0;
for k=1:length(current0)
    if abs(voltage-current0(k,1)) < 0.25
        currFactor = 1.0/current0(k,2);
        break
    end
end
beta = 1./logspace(expTempInitial, expTempFinal, nAttempts);

inFile = sprintf('%s/%s.dat', inDir, inPrefix);
outFile = sprintf('%s/%s_comb.dat', outDir, inPrefix);
data = dlmread(inFile, ' ');
x = data(:,1);
y = data(:,2);
n = length(x);

disp(sprintf('\nTrace monte'))
% Define the initial partitioning.
%bestErr = std(y)/sqrt(n);
bestErr = 1e20;
bestPart = [ones(nBins,1); n];

lastErr = bestErr;
partition = ceil((rand(nBins-1,1) + 1e-9)*n);
for j=1:nAttempts
    % Randomly chose partitions and sort them.
    shift = round(step*(rand(nBins-1,1) - 0.5));
    part = [1; partition+shift; n];
    part = sort(part);
    
    % Compute the length of the partitions.
    num = diff(part);
    
    % Compute the standard error in the bins from the partitions.
    err = zeros(nBins,1);
    nActive = 0;
    for j=1:nBins   
        if num(j) <= 1, continue, end
        
        %err(j) = std(y(part(j):part(j+1)))/sqrt(num(j));
        err(j) = std(y(part(j):part(j+1)))/sqrt(num(j));
        nActive = nActive + 1;
    end
    currErr = sum(err.^2);
    
    % Metropolis
    dErr = (currErr - lastErr)*beta(j);
    if nActive == nBins && (currErr < lastErr || rand < exp(-dErr))
        partition = partition + shift;
        lastErr = currErr;
                
        if currErr < bestErr
            bestErr = currErr;
            bestPart = part;
            disp(sprintf('Best error: %f', currErr))
        end
    end
end

% Compute values for the best partition.
part = bestPart;
num = diff(part);
meanX = zeros(nBins,1);
meanY = zeros(nBins,1);
seY = zeros(nBins,1);
for j=1:nBins
    meanX(j) = mean(x(part(j):part(j+1)));
    meanY(j) = mean(y(part(j):part(j+1)));
    
    if num(j) <= 1
        seY(j) = 0;
    else
        seY(j) = std(y(part(j):part(j+1)))/sqrt(num(j));
    end
end

errorbar(meanX, meanY*currFactor, seY*currFactor)

dlmwrite(outFile, [meanX meanY*currFactor seY*currFactor], ' ');




