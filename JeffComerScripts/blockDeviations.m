% blockDeviations.m
function [block blockDev] = blockDeviations(t, x, minBlockSize, n)

blockSize = logspace(log10(minBlockSize), log10(0.5*length(x)), n);
blockSize = unique(round(blockSize));
blocks = length(blockSize);
block = zeros(blocks,1);
blockDev = zeros(blocks,1);

for j=1:length(blockSize)
  block(j) = blockSize(j);
  [tavg xavg xerr] = blockAverage(t, x, blockSize(j));
  blockDev(j) = std(xavg)/sqrt(length(xavg));
end
