clear all;

% Parameters:
dir = 0;
lx = 11.0; % in nm (the input data is converted to nm also)
ly = 11.0; % in nm
lz = 25.0; % in nm
rangeV = [-3 1.2]; % voltage range (V)
numLevelsV = 40; % number of contours to draw
colorScale = [0.1 0.1 0.4; 0.2 0.2 0.6; 0.3 0.3 0.8; 0.4 0.4 1; 1 1 1; 1 0 0; 0.8 0 0; 0.6 0 0; 0.4 0 0]; 
% Input:
prefix = 'cross_';
nameList = {'dopc_q2c'};
% Output:
outName = 'plot_';

% Tick marks
tickX = [-8 -6 -4 -2 0 2 4 6 8];
tickY = [-8 -6 -4 -2 0 2 4 6 8];
fancyTicks = 0; % Use Jeff's fancy tickmarks?

% Controls the portion displayed
lFactor = 1.0;

for pick=1:length(nameList)
    % Input:
    name = nameList{pick};
    crossFile = sprintf('%s%s.dat', prefix, name);
    crossFileX = sprintf('%s%s.dat.rx', prefix, name);
    crossFileY = sprintf('%s%s.dat.ry', prefix, name);
    crossFileZ = sprintf('%s%s.dat.rz', prefix, name);
    % Output:
    crossPlot = sprintf('%s%s.pdf', outName, name);
    levV = linspace(rangeV(1),rangeV(2),numLevelsV);
    
    % Make the color map.
    mapV = makeColorMap(colorScale, length(levV)+1);
    
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
    
    % Plot the contour map.
    h = figure(pick);
    set(h, 'Position', [10 10 500 800], 'Resize', 'off');
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
    
    % Add Jeff's fancy tickmarks.
    if fancyTicks
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
    
    % Generate a pdf of the plot.
    print('-dpdf', crossPlot)
end
