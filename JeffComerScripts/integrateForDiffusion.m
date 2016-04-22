clear all;

prefix = 'quick';
kT = 0.5862292;
velConvFactor = 15*15*1e6; % From A^2/fs to A^2/ns
forceConvFactor = 1e6; % From A^2/fs to A^2/ns
resolution = 60;

ghost = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the velocity autocorrelation function.
% Load the data.
data = dlmread(strcat(prefix,'.vcorr'), ' ');
velTime = data(:,1);
velCorr = data(:,2);

% Add ghost points.
velDt = velTime(2)-velTime(1);
for j=1:ghost
    velTime = [velTime(1)-velDt; velTime; velTime(end)+velDt];
    velCorr = [velCorr(1); velCorr; velCorr(end)];
end

% Form the high resolution time.
velN = length(velTime);
velHighTime = linspace(velTime(1),velTime(end)-velDt/resolution,resolution*velN);

% Compute the splines.
velHermite = pchip(velTime, velCorr, velHighTime);
velSpline = spline(velTime, velCorr, velHighTime);

% Remove the ghost points.
cut = resolution*ghost;
velTime = velTime(1+ghost:end-ghost);
velCorr = velCorr(1+ghost:end-ghost);
velHighTime = velHighTime(1+cut:end-cut);
velHermite = velHermite(1+cut:end-cut);
velSpline = velSpline(1+cut:end-cut);

% Compute the diffusivity.
velDiffRaw = velConvFactor*trapz(velTime, velCorr)/3;
velDiffHermite = velConvFactor*trapz(velHighTime, velHermite)/3;
velDiffSpline = velConvFactor*trapz(velHighTime, velSpline)/3;
velIntRaw = cumtrapz(velTime, velCorr);
velIntHermite = cumtrapz(velHighTime, velHermite);
velIntSpline = cumtrapz(velHighTime, velSpline);

% Write the results.
fprintf('\n');
fprintf('velDiffRaw: %g\n', velDiffRaw);
fprintf('velDiffHermite: %g\n', velDiffHermite);
fprintf('velDiffSpline: %g\n', velDiffSpline);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the force autocorrelation function.
% Load the data.
data = dlmread(strcat(prefix,'.fcorr'), ' ');
forceTime = data(:,1);
forceCorr = data(:,2);

% Add ghost points.
forceDt = forceTime(2)-forceTime(1);
for j=1:ghost
    forceTime = [forceTime(1)-forceDt; forceTime; forceTime(end)+forceDt];
    forceCorr = [forceCorr(1); forceCorr; forceCorr(end)];
end

% Form the high resolution time.
forceN = length(forceTime);
forceHighTime = linspace(forceTime(1),forceTime(end)-forceDt/resolution,resolution*forceN);

% Compute the splines.
forceHermite = pchip(forceTime, forceCorr, forceHighTime);
forceSpline = spline(forceTime, forceCorr, forceHighTime);

% Remove the ghost points.
cut = resolution*ghost;
forceTime = forceTime(1+ghost:end-ghost);
forceCorr = forceCorr(1+ghost:end-ghost);
forceHighTime = forceHighTime(1+cut:end-cut);
forceHermite = forceHermite(1+cut:end-cut);
forceSpline = forceSpline(1+cut:end-cut);

% Compute the diffusivity.
forceDiffRaw = forceConvFactor*kT^2/(sum(forceCorr)*forceDt);
%forceDiffRaw = forceConvFactor*kT^2/trapz(forceTime, forceCorr);
forceDiffHermite = forceConvFactor*kT^2/trapz(forceHighTime, forceHermite);
forceDiffSpline = forceConvFactor*kT^2/trapz(forceHighTime, forceSpline);
forceIntRaw = cumtrapz(forceTime, forceCorr);
forceIntHermite = cumtrapz(forceHighTime, forceHermite);
forceIntSpline = cumtrapz(forceHighTime, forceSpline);

fprintf('\n');
fprintf('forceDiffRaw: %g\n', forceDiffRaw);
fprintf('forceDiffHermite: %g\n', forceDiffHermite);
fprintf('forceDiffSpline: %g\n', forceDiffSpline);

figure(1)
plot(velTime, velCorr, 'ko', velHighTime, velHermite, 'r-', velHighTime, velSpline, 'b-')
axis([0 2000 -inf inf])
xlabel('time (fs)', 'FontSize', 20)
ylabel('velocity autocorr. (A^2/fs^2)', 'FontSize', 20)
%saveas(gcf,  strcat(prefix,'.vel.png'), 'png');

figure(2)
plot(forceTime, forceCorr, 'ko', forceHighTime, forceHermite, 'r-', forceHighTime, forceSpline, 'b-')
axis([0 2000 -inf inf])
xlabel('time (fs)', 'FontSize', 20)
ylabel('force autocorr. (kT^2/A^2)', 'FontSize', 20)
%saveas(gcf,  strcat(prefix,'.force.png'), 'png');

figure(3)
plot(velTime, velIntRaw, 'ko', velHighTime, velIntHermite, 'r-', velHighTime, velIntSpline, 'b-')
axis([0 2000 -inf inf])
xlabel('time (fs)', 'FontSize', 20)
ylabel('integrated vel. autocorr. (A/fs)^2*fs', 'FontSize', 16)
legend('trapezoid', 'Hermite', 'spline');
saveas(gcf,  strcat(prefix,'.vel.png'), 'png');

figure(4)
plot(forceTime, forceIntRaw, 'ko', forceHighTime, forceIntHermite, 'r-', forceHighTime, forceIntSpline, 'b-')
axis([0 2000 -inf inf])
xlabel('time (fs)', 'FontSize', 20)
ylabel('integrated for. autocorr. (kcal/(mol A))^2 fs', 'FontSize', 16)
legend('trapezoid', 'Hermite', 'spline');
saveas(gcf,  strcat(prefix,'.force.png'), 'png');
