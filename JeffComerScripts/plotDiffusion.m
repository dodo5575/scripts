clear all;
% Author: Jeff Comer <jcomer@illinois.edu>

% Parameters:
corrPrefix = 'correlation/diffuse_no_1M_pot';
%corrPrefix = 'correlation/basepair_no_chl';
minInstances = 200;
% More parameters.
corrSuffix = 'dat.vcorr';
logSuffix = 'dat.log';
corrFormat = '%s.%d.%s';
cylFile = 'cyl2.dat';
rangeV = [100 300];
nLevels = 8;
tick = 0:5:35;
velConvFactor = 1e6; % From A^2/fs to A^2/ns

% Load the cylindrical grid's geometry.
data = dlmread(cylFile, ' ');
s0 = data(1,1);
ds = data(2,1);
ns = data(3,1);
z0 = data(1,2);
dz = data(2,2);
nz = data(3,2);
nGrid = ns*nz;

% Load the velocity autocorrelation function if it exists.
diffuse = zeros(nGrid,1);
nValid = 0;
totalInst = 0;
usedInst = 0;
for j=0:nGrid-1
    inst = 0;
    logName = sprintf(corrFormat, corrPrefix, j, logSuffix);
    if exist(logName, 'file')
        inst = readWeight(logName);
        totalInst = totalInst + inst;
        if inst < minInstances
            %fprintf('node %d: %d instances\n', j, inst);
        else
            usedInst = usedInst + inst;
        end
    else
        %fprintf('node %d: does not exist\n', j);
    end
    
    % Find the velocity correlation function file corresponding to this
    % node.
    % Those that don't exist or those that do not represent averages
    % over a sufficient number of inst.
    fileName = sprintf(corrFormat, corrPrefix, j, corrSuffix);
    if exist(fileName, 'file') && inst >= minInstances
        % Read the file.
        data = dlmread(fileName, ' ');
        %fprintf('Read %s.\n', fileName);
        
        % Integrate the autocorrelation function.
        diffuse(j+1) = velConvFactor*computeDiffusivity(data(:,1), data(:,2));
        nValid = nValid + 1;
    end
end
fprintf('%d of %d nodes are valid\n', nValid, nGrid);
fprintf('%d of %d instances are used\n', usedInst, totalInst);

% Store the positions of each node.
s = zeros(nGrid,1);
z = zeros(nGrid,1);
v = zeros(nGrid,1);
j = 1;
for is=0:ns-1
  for iz=0:nz-1
    s(j) = s0 + is*ds;
    z(j) = z0 + iz*dz;
    v(j) = diffuse(j);
    %fprintf('%g %g\n', s, z);
    j = j + 1;
  end
end

%mat = regexp(name, '[0123456789\.]+V', 'match')
% Prepare the plot.
levV = linspace(rangeV(1),rangeV(2),nLevels);
mapV = makeColorMap([0 0 0; 0.15 0.15 1; 0.6 1 0.6; 0.9 1 0.9; 1 1 1], length(levV)+1);

% Make the section.
secX = zeros(ns,nz);
secY = zeros(ns,nz);
secV = zeros(ns,nz);
for is=0:ns-1
  for iz=0:nz-1
    j = 1 + iz + is*nz;
    secX(is+1, iz+1) = s(j);
    secY(is+1, iz+1) = z(j);
    secV(is+1, iz+1) = v(j);
  end
end

figure(1)
pcolor(secX, secY, secV)
%contourf(secX, secY, secV, levV)
colormap(mapV)
xlabel('s (A)', 'FontSize', 16)
ylabel('z (A)', 'FontSize', 16)
%axis vis3d
%shading('interp')
axis square
caxis(rangeV)
colorbar
set(gcf, 'color', 'white'); % sets the color to white
set(gca, 'TickDir', 'in', 'XTick', tick, 'YTick', tick);
set(gca, 'FontSize', 16);
tickSize = 0.03*(secX(end)-secX(1));
tickSize1 = 0.035*(secX(end)-secX(1));
for j=1:length(tick)
    h = line([tick(j) tick(j)], [secY(1) secY(1)+tickSize1]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([tick(j) tick(j)], [secY(1) secY(1)+tickSize]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
    
    h = line([tick(j) tick(j)], [secY(end)-tickSize1 secY(end)]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([tick(j) tick(j)], [secY(end)-tickSize secY(end)]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
    
    h = line([secX(1) secX(1)+tickSize1], [tick(j) tick(j)]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([secX(1) secX(1)+tickSize], [tick(j) tick(j)]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
    
    h = line([secX(end)-tickSize1 secX(end)], [tick(j) tick(j)]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([secX(end)-tickSize secX(end)], [tick(j) tick(j)]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
end
%print('-dpdf', depthPlot)
