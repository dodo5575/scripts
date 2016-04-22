clear all;
inFile = 'pos_trap_0V0-0.dat';

kB = 1.380658e-23; % J/K
temp = 295; % K

data = dlmread(inFile, ' ');
z = data(:,2);

zVar = var(z)
zMean = mean(z)
zSqDev = (z-zMean).^2;
zVar1 = mean(zSqDev)
zVarSe = std((z-zMean).^2)/length(zSqDev)

springMean = kB*temp/(zVar*1e-18)
springErr = kB*temp/(zVar^2*1e-18)*zVarSe

[binN binZ] = hist(z,20);
binSize = binZ(2)-binZ(1);
prob = binN./sum(binSize.*binN);

zSort = sortrows(z);
estProb = 1/sqrt(2*pi*zVar)*exp(-(zSort-zMean).^2/(2*zVar));

figure(1)
gh = plot(binZ, prob,'bo-', zSort, estProb, 'k-');
set(gh, 'MarkerSize', 8);
set(gh, 'LineStyle', '-');
    
xlabel('z (nm)')
ylabel('number')



