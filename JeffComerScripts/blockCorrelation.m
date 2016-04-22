clear all;

name = 'spacingOne1';
%name = 'spacing5';

inFile = sprintf('%s.log.pos', name);

data = dlmread(inFile, ' ');
t = data(:,1);
z = data(:,2);
%z = rand(size(t));
dt = t(2) - t(1);

[block blockDev] = blockDeviations(t, z, 20, 500);
t = block*dt;

plot(t, blockDev, 'k-');

outFile = sprintf('blockStd_%s.dat', name);
dlmwrite(outFile, [t' blockDev'], ' ');
