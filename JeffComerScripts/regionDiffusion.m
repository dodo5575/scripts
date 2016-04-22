clear all;
% Author: Jeff Comer <jcomer@illinois.edu>

nameList = {'diffuse_no_1M' 'basepair_no'};
ionList = {'pot' 'chl'};
%ionList = {'pot'};
inDir = 'correlation_short';
outDir = 'region_integral';

% Parameters:
velConvFactor = 1e6; % From A^2/fs to A^2/ns
nSteps = 60; % Number of points in each instance
% Input:
diffusionGrid = 'diffGrid.dx';
cylFile = 'cyl2.dat';
corrSuffix = 'dat.vcorr';
logSuffix = 'dat.log';
corrFormat = '%s.%d.%s';
% Region parameters.
sIn = 3;
zIn = 15;
zTop = 21;
sTop = 11;
zIntercept = 13;
% Output:
outPrefix = 'region';

% Load the cylindrical grid's geometry.
data = dlmread(cylFile, ' ');
s0 = data(1,1);
ds = data(2,1);
ns = data(3,1);
z0 = data(1,2);
dz = data(2,2);
nz = data(3,2);
nGrid = ns*nz;
plotColor = [0 0 0; 0.9 0 0; 0.8 0.4 0; 0.5 0.5 0; 0 0.8 0; 0 0.7 0.7; 0 0 0.9; 0.6 0 0.6; 0.6 0.6 0.6];

figure(1)
clf
hold on
fprintf('\n****%s\n', datestr(now));

count = 0;
plotLegend = {};
for ionInd=1:length(ionList)
    for nameInd=1:length(nameList)
        name = nameList{nameInd};
        ion = ionList{ionInd};
        corrPrefix = sprintf('%s/%s_%s', inDir, name, ion);
        
        % Load the velocity autocorrelation function if it exists.
        corrVel = zeros(nSteps,nGrid);
        corrTime = zeros(nSteps,nGrid);
        weight = zeros(nGrid,1);
        totalInst = 0;
        for j=0:nGrid-1
            inst = 0;
            % Find the velocity correlation function file corresponding to this
            % node.
            logName = sprintf(corrFormat, corrPrefix, j, logSuffix);
            fileName = sprintf(corrFormat, corrPrefix, j, corrSuffix);
            
            if exist(logName, 'file') && exist(fileName, 'file')
                inst = readWeight(logName);
                totalInst = totalInst + inst;
                weight(j+1) = inst;
                
                % Read the file.
                data = dlmread(fileName, ' ');
                %fprintf('Read %s.\n', fileName)
                corrTime(:,j+1) = data(:,1);
                corrVel(:,j+1) = data(:,2);
            end
        end
        fprintf('\nSystem %s_%s.\n', name, ion);
        fprintf('Read %d instances.\n', totalInst);
        
        % Store the positions of each node.
        % Find the region that each pertain to.
        s = zeros(nGrid,1);
        z = zeros(nGrid,1);
        reg = zeros(nGrid,1);
        j = 1;
        for is=0:ns-1
            for iz=0:nz-1
                s(j) = s0 + is*ds;
                z(j) = z0 + iz*dz;
                
                % Find which region we are in.
                where = 0;
                if z(j) <= zIn && s(j) <= sIn
                    where = 2;
                elseif z(j) <= zIn && s(j) > sIn
                    where = 3;
                elseif z(j) <= zTop && s(j) > sTop
                    where = 3;
                elseif z(j) > zTop
                    where = 1;
                elseif z(j) > zTop && s(j) <= sIn
                    where = 1;
                end
                
                % Handle the region for the corner of the pore.
                if where == 0
                    if z(j) > 1.0*s(j) + 13
                        where = 1;
                    else
                        where = 3;
                    end
                end
                
                reg(j) = where;
                %fprintf('%g %g\n', s, z);
                j = j + 1;
            end
        end
        
        % Make the average correlation function for each region.
        regionName = {'out' 'pore' 'surf'};
        for k=1:3
            t = zeros(nSteps,1);
            c = zeros(nSteps,1);
            w = 0;
            for j=1:nGrid
                if reg(j) == k
                    t = corrTime(:,j);
                    c = c + corrVel(:,j).*weight(j);
                    w = w + weight(j);
                end
            end
            
            % Normalize.
            c = c/w;
            diffuse = velConvFactor*computeDiffusivity(t, c);
            [intT intDiff] = integrateDiffusivity(t, c);
            intT = intT*velConvFactor;
            intDiff = intDiff*velConvFactor;
            
            fprintf('%s(%d): %g\n', regionName{k}, w, diffuse);
            
            % Plot it.
            count = count + 1;
            color0 = plotColor(mod(count-1,length(plotColor))+1,:);
            gh = plot(intT,intDiff);
            set(gh, 'Color', color0);
            set(gh, 'Marker', 'None');
            set(gh, 'MarkerSize', 7);
            set(gh, 'LineStyle', '-');
            plotLegend{count} = sprintf('%s_%s_%s', regionName{k}, name, ion);
            
            % Write it.
            outFile = sprintf('%s/%s_%s_%s_%s.int', outDir, outPrefix, regionName{k}, name, ion);
            dlmwrite(outFile, [intT intDiff], ' ');
            %fprintf('Wrote "%s".\n', outFile);
        end    
    end
end

hold off
h = legend(plotLegend);
set(h, 'Interpreter', 'none');
