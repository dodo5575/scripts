clear all;

dataNum = 0;
dataRun = 'realBond';
dataName = sprintf('/projects/jcomer/pmf/vacuum/analysis/data/three_%s_chl-chl', dataRun);
dataR0 = 0.3*(dataNum + 1);

kT = 0.5862292;
spring = 6.4;
n = 100000;

r0 = dataR0;
rc = 0.5*r0 + 0.5*sqrt(r0^2 + 8*kT/spring);
rMax = 3*rc;

r = linspace(0.0, rMax, 600)';
dr = r(2)-r(1);

% Show the normalized distribution.
analProb = 4*pi*r.^2.*exp(-0.5*spring*(r-r0).^2/kT);
norma = trapz(r,analProb);
analProb = analProb/norma;

rMean = trapz(r, analProb.*r)
diffMean = rMean - r0

% Generate the random distribution.
samp = bondDistribRandom(n, r0, spring, kT);
[sampProb sampR] = hist(samp, 60);
norma = trapz(sampR,sampProb);
sampProb = sampProb/norma;

% Extract the distribution from the file.
inFile = sprintf('%s%d.dat', dataName, dataNum);
data = dlmread(inFile, ' ');
[dataProb dataR] = hist(data(:,4), 30);
norma = trapz(dataR,dataProb);
dataProb = dataProb/norma;

figure(1)
clf
%plot(sampR, sampProb, 'r--', dataR, dataProb, 'bo-')
plot(r,analProb,'k--', sampR, sampProb, 'r.-', dataR, dataProb, 'bo-')
%plot(r,analProb,'k--', dataR, dataProb, 'b-')

