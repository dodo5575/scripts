clear all;

inPrefix = 'curr_vs_com_double';
binSize = 8;

inFile = sprintf('%s.dat', inPrefix);
outFile = sprintf('%s_bin.dat', inPrefix);
data = dlmread(inFile, ' ');
x = data(:,1);
y = data(:,2);
n = length(x);

% Sort the data and average in bins.
dataSort = sortrows([x y]);
nBins = length(1:binSize:n);
meanX = zeros(nBins,1);
seX = zeros(size(meanX));
meanY = zeros(size(meanX));
seY = zeros(size(meanX));

for j=1:nBins
    i0 = 1 + (j-1)*binSize;
    i1 = j*binSize;
    
    % Handle the end of the list.
    if i1 > n
        i1 = n;
    end
    
    nSamples = length(i0:i1);
    
    meanX(j) = mean(dataSort(i0:i1,1));
    seX(j) = std(dataSort(i0:i1,1)) / sqrt(nSamples);
    %seX(j) = std(dataSort(i0:i1,1));
    meanY(j) = mean(dataSort(i0:i1,2));
    seY(j) = std(dataSort(i0:i1,2)) / sqrt(nSamples);
    %seY(j) = std(dataSort(i0:i1,2));
end

figure(1)
errorbar(meanX, meanY, seY)

dlmwrite(outFile, [meanX meanY seY], ' ');





