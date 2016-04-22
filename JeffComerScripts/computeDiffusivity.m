function diffuse = computeDiffusivity(velTime, velCorr) 
% Compute the diffusivity given a velocity autocorrelation function.
% Author: Jeff Comer <jcomer2@illinois.edu>
ghost = 4;
resolution = 60;

% Add ghost points.
velDt = velTime(2)-velTime(1);
for j=1:ghost
    velTime = [velTime(1)-velDt; velTime; velTime(end)+velDt];
    velCorr = [velCorr(1); velCorr; velCorr(end)];
end

% Form the high resolution time.
velN = length(velTime);
velHighTime = linspace(velTime(1),velTime(end)-velDt/resolution, resolution*velN);

% Compute the spline.
velHermite = pchip(velTime, velCorr, velHighTime);

% Remove the ghost points.
cut = resolution*ghost;
velTime = velTime(1+ghost:end-ghost);
velCorr = velCorr(1+ghost:end-ghost);
velHighTime = velHighTime(1+cut:end-cut);
velHermite = velHermite(1+cut:end-cut);

% Compute the diffusivity.
%velDiffRaw = trapz(velTime, velCorr)/3;
velDiffHermite = trapz(velHighTime, velHermite)/3;
%velIntRaw = cumtrapz(velTime, velCorr);
%velIntHermite = cumtrapz(velHighTime, velHermite);

diffuse = velDiffHermite;
