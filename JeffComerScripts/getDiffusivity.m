clear all;

inDir = '.';
fileGlob = '*_dopc*.msd';
cutTime = 1.0;
dimensions = 2;

plotColor = [0 0 0; 0.9 0 0; 0.8 0.4 0; 0.5 0.5 0; 0 0.8 0; 0 0.7 0.7; 0 0 0.9; 0.6 0 0.6; 0.6 0.6 0.6];

figure(1)
clf
hold on

glob = sprintf('%s/%s', inDir, fileGlob);
fileList = dir(glob);
nFiles = length(fileList);
fprintf('Found %d files.\n', nFiles);

for j=1:nFiles
    fileName = sprintf('%s/%s', inDir, fileList(j).name);    
    data = dlmread(fileName, ' ');
    
    dt = data(2,1) - data(1,1);
    cutIndex = floor(cutTime/dt);
    
    if length(data(:,1)) < cutIndex
        continue
    end
    
    % Get the data.
    t = data(cutIndex:end,1);
    n = length(t);
    msd = data(cutIndex:end,2);
    
    % Compute the linear regression.
    [a b da db] = linearRegression(t, msd);
    theo = b*t + a;
    
    % Compute the diffusivity.
    diffusion = b/(2*dimensions);
    diffusionErr = db/(2*dimensions);
    
    % Plot it.
    color0 = plotColor(mod(j-1,length(plotColor))+1,:);
    gh = plot(t, msd);
    set(gh, 'Color', color0);
     set(gh, 'Marker', 'o');
    set(gh, 'MarkerSize', 7);
    set(gh, 'LineStyle', 'None');
   
    gh = plot(t, theo);
    set(gh, 'Color', color0);
    set(gh, 'LineStyle', '-');
    
    % Say the result.
    dispText = sprintf('slope = %.3f, D = %.3f +/- %.3f', b, diffusion, diffusionErr);
    %text(0.1, 4+0.5*j, dispText);
    fprintf('diffusivity: %s %.3f +/- %.3f\n', fileName, diffusion, diffusionErr);
end

xlabel('t (ns)');
ylabel('<x^2 + y^2> A^2');
hold off

