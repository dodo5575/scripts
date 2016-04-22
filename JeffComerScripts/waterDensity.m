clear all;

sysList = {'b-dna' 's-dna' 'open'};
dcdFreq = 5000*20;
hexVec1 = [129.115 0 0];
hexVec2 = [64.5575 111.817 0];
areaTotal = norm(cross(hexVec1,hexVec2));
rad = 30;
z0 = 70;
z1 = 125;
area = areaTotal - pi*rad*rad;
vol = 2*(z1-z0)*area;
cutTime = 3.5;
unifiedMass = 1.660538782e-27;

cutInd = ceil(cutTime/dcdFreq);
fprintf('sys\t\tdensity(kg/m^3)\t\terr(kg/m^3)\n');

for s=1:length(sysList)
    inFile = sprintf('numWater_pore6_%s_500mV_far.dat', sysList{s});
    data = dlmread(inFile, ' ');
    numWater = data(cutInd:end,2);
    
    nw = mean(numWater);
    ew = std(numWater)/sqrt(length(numWater));
    density = (unifiedMass*1e30)*18.0106*nw/vol;
    densityErr = (unifiedMass*1e30)*18.0106*ew/vol;
    
    %fprintf(' %g +/- %g\n', nw, ew);
    fprintf('%s\t\t%g\t\t%g\n', sysList{s}, density, densityErr);
end
