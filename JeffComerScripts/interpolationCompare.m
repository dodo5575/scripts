clear all;

stride = 20;
corrFuncName = 'diffusion_chl_farCorrV.dat';
ghost = 4;
force = 0;
kT = 0.5862292;

data = dlmread(corrFuncName, ' ');
rawT = data(:,1);
rawF = data(:,2);
rawN = length(rawT);

dataT = rawT(1:stride:end);
dataF = rawF(1:stride:end);
dataN = length(dataT);

highT = rawT;

% Add ghost points.
dt = dataT(2)-dataT(1);
for j=1:ghost
    dataT = [dataT(1)-dt; dataT; dataT(end)+dt];
    dataF = [dataF(1); dataF; dataF(end)];
end

% Add the high res ghost points.
dt = rawT(2)-rawT(1);
for j=1:stride*ghost
    highT = [highT(1)-dt; highT; highT(end)+dt];
end

% Compute the splines.
hermiteF = pchip(dataT, dataF, highT);
splineF = spline(dataT, dataF, highT);
cut = stride*ghost;

% Remove the ghost points.
dataT = dataT(1+ghost:end-ghost);
dataF = dataF(1+ghost:end-ghost);
highT = highT(1+cut:end-cut);
hermiteF = hermiteF(1+cut:end-cut);
splineF = splineF(1+cut:end-cut);

if force
    rawArea = 1e6*kT^2/trapz(rawT, rawF);
    dataArea = 1e6*kT^2/trapz(dataT, dataF);
    splineArea = 1e6*kT^2/trapz(highT, splineF);
    hermiteArea = 1e6*kT^2/trapz(highT, hermiteF);
else
    rawArea = 1e6*trapz(rawT, rawF)/3;
    dataArea = 1e6*trapz(dataT, dataF)/3;
    splineArea = 1e6*trapz(highT, splineF)/3;
    hermiteArea = 1e6*trapz(highT, hermiteF)/3;
end

fprintf('\n');
fprintf('rawArea: %g\n', rawArea);
fprintf('dataArea: %g\n', dataArea);
fprintf('splineArea: %g\n', splineArea);
fprintf('hermiteArea: %g\n', hermiteArea);

figure(1)
plot(rawT, rawF, 'k--', dataT, dataF, 'ro', highT, hermiteF, 'b-', highT, splineF, 'g--');
legend('raw', 'data', 'hermite', 'spline')
%axis([ inf -inf inf])
