clear all;

dcdFreq = 5000; % in fs
histBins = 16;

cutTime = 6; % in ns
name = 'trap2.0';
inFile = 'pos_cen/pos_trap2.0_0V.dat';
corrSamples = 100;

outHist = sprintf('hist_%s_0Va.dat', name);
outMod = sprintf('model_%s_0Va.dat', name);

boltz = 1.3806504e-23;
temp = 295;
kT = boltz*temp; % in J

cut = floor(cutTime*1e6/dcdFreq);
data = dlmread(inFile, ' ');
z = data(cut:end,2);

meanZ = mean(z)
stdZ = std(z)
z = z - meanZ;

[histN histZ] = hist(z,histBins);
histN = histN';
histZ = histZ';
dz = histZ(2)-histZ(1);
histProb = histN/(length(z)*dz);

correl = xcorr(z);
mid = floor(length(correl)/2);
correl = correl(mid:end);
corrLog = log(correl);
corrN = length(correl);
corrSamp = (1:corrN)';
dlmwrite(sprintf('logcorr_%s.dat', name), corrLog, ' ');

modZ = linspace(histZ(1),histZ(end),400)';
modProb = 1/(stdZ*sqrt(2*pi))*exp(-0.5*(modZ/stdZ).^2);
modEff = 1e18*kT/(stdZ*stdZ); %/home/jcomer/projects/periodic_trap/analysis in nN/nm
modN = mean(histN)/mean(modProb)*modProb;
expectProb = 1/(stdZ*sqrt(2*pi))*exp(-0.5*(histZ/stdZ).^2);

% Compute the correlation time.
delta = corrN*sum(corrSamp.*corrSamp) - sum(corrSamp)^2; 
slope = (sum(corrSamp.*corrSamp)*sum(corrLog) - sum(corrSamp)*sum(corrSamp.*corrLog))/delta;
offset = (corrN*sum(corrSamp.*corrLog) - sum(corrSamp)*sum(corrLog))/delta;
corrFit = slope*corrSamp + offset;
%corrSamples = floor(1/slope)

%chiSq = sum((expectProb-histProb).^2./expectProb)*length(expectProb)/2;
indepSamples = length(z)/corrSamples;
modErr = modEff/sqrt(2*(indepSamples-1));

%plot(histZ,histProb,'bo-',modZ,modProb,'k--')
figure(1)
plot(histZ,histN,'bo-',modZ,modN,'k--')

figure(2)
%plot(corrSamp,corrLog,'ko',corrSamp,corrFit,'b--')
plot(corrSamp,correl,'ko')

disp(sprintf('Efficiency: %f +/- %f nN/nm', modEff, modErr));
%disp(sprintf('Reduced chi squared: %f', chiSq));

dlmwrite(outHist,[histZ histProb], ' ');
dlmwrite(outMod,[modZ modProb], ' ');

