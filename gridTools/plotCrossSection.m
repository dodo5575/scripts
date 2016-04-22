clear all;

% Parameters:
dir = 0;
prefix = 'cross/sym_conc';
outName = 'plot';
%nameList = {'K_open' 'K_b-dna' 'K_s-dna' 'Cl_open' 'Cl_b-dna' 'Cl_s-dna'};
%nameList = {'K_b-dna' 'Cl_b-dna'};
nameList = {'K_s-dna' 'Cl_s-dna'};
%nameList = {'K_open' 'Cl_open'};

rangeV = [-0.02 0.6];
%levV = [];
lFactor = 0.8;
lx = 12.0;
ly = 12.0;
lz = 25.0;
tickX = [-4 -2 0 2 4];
tickY = [-8 -4 0 4 8];
plotType = 'png';

for pick=1:length(nameList)
    % Input:
    name = nameList{pick};
    crossFile = sprintf('%s%s.dat', prefix, name);
    crossFileX = sprintf('%s%s.dat.rx', prefix, name);
    crossFileY = sprintf('%s%s.dat.ry', prefix, name);
    crossFileZ = sprintf('%s%s.dat.rz', prefix, name);
    % Output:
    crossPlot = sprintf('%s%s.%s', outName, name, plotType);
    
    levV = linspace(rangeV(1),rangeV(2),12);
    
    %mapV = makeColorMap([0.4 0 0; 0.6 0 0; 0.8 0 0; 1 0 0; 1 1 1; 1 1 1; 0 0 1; 0 0 0.8; 0 0 0.6; 0 0 0.4], 3*length(levV));
    
    mapV = makeColorMap([0 0 0; 0.3 0.3 1; 0.6 1 0.6; 0.9 1 0.9; 1 1 1], length(levV)+1);
    
    
    %mapV = makeColorMap([0 0 0; 0.4 0.34 0.28; 0.8 0.68 0.55; 1 1 0.85; 1 1 1], length(levV)+1);
    
    % Read the data.
    v = dlmread(crossFile, ' ');
    x = dlmread(crossFileX, ' ');
    y = dlmread(crossFileY, ' ');
    z = dlmread(crossFileZ, ' ');
    
    % Convert to nanometers from angstroms.
    z = 0.1*z;
    x = 0.1*x;
    y = 0.1*y;
    
    if dir == 0
        cenX = y;
        cenY = z;
    elseif dir == 1
        cenX = x;
        cenY = z;
    else
        cenX = x;
        cenY = y;
    end
    cenV = v;
    
    h = figure(pick);
    set(h, 'Position', [10 10 420 800], 'Resize', 'off');
    %set(h, 'Position', [10 10 500 800]);
    %pcolor(cenX, cenY, cenV)
    %shading('interp')
    contourf(cenX, cenY, cenV, levV)
    colormap(mapV)
    xlabel('x (nm)', 'FontSize', 20)
    ylabel('y (nm)', 'FontSize', 20)
    axis([-0.5*lFactor*lx 0.5*lFactor*lx -0.5*lFactor*ly 0.5*lFactor*lz])
    axis equal
    %axis square
    caxis(rangeV)
    colorbar
    set(gcf, 'color', 'white'); % sets the color to white
    set(gca, 'TickDir', 'out', 'XTick', tickX, 'YTick', tickY, 'TickLength', [0.02 0.02]);
    set(gca, 'FontSize', 20);
    if 0
        tickSize = 0.03*(cenX(end)-cenX(1));
        tickSize1 = 0.035*(cenX(end)-cenX(1));
        for j=1:length(tickX)
            h = line([cenX(1) cenX(1)+tickSize1], [tickX(j) tickX(j)]);
            set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
            h = line([cenX(1) cenX(1)+tickSize], [tickX(j) tickX(j)]);
            set(h,'LineWidth', 2, 'Color', [0 0 0]);
            
            h = line([cenX(end)-tickSize1 cenX(end)], [tickX(j) tickX(j)]);
            set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
            h = line([cenX(end)-tickSize cenX(end)], [tickX(j) tickX(j)]);
            set(h,'LineWidth', 2, 'Color', [0 0 0]);
        end
        for j=1:length(tickY)
            h = line([tickY(j) tickY(j)], [cenY(1) cenY(1)+tickSize1]);
            set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
            h = line([tickY(j) tickY(j)], [cenY(1) cenY(1)+tickSize]);
            set(h,'LineWidth', 2, 'Color', [0 0 0]);
            
            h = line([tickY(j) tickY(j)], [cenY(end)-tickSize1 cenY(end)]);
            set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
            h = line([tickY(j) tickY(j)], [cenY(end)-tickSize cenY(end)]);
            set(h,'LineWidth', 2, 'Color', [0 0 0]);
        end
    end
    
    saveas(gcf, crossPlot, plotType);
end
