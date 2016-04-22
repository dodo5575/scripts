clear all;

% Parameters:
outName = 'first_';
dir = 2;
run = 'crossZ_conc_POT_';
nameList = {'ade_conc0.1' 'ade_brown_diel176_conc0.1' 'ade_conc1.0' 'ade_brown_diel176_conc1.0'};
rangeV = [0 8];
levN = 20;

lFactor = 1.0;
lx = 5;
ly = 5;
tickX = [-1 0 1];
tickY = [-1 0 1];

for pick=1:length(nameList)
    % Input:
    name = nameList{pick};
    crossFile = sprintf('%s%s.dat', run, name);
    crossFileX = sprintf('%s%s.dat.rx', run, name);
    crossFileY = sprintf('%s%s.dat.ry', run, name);
    crossFileZ = sprintf('%s%s.dat.rz', run, name);
    % Output:
    crossPlot = sprintf('%s%s', outName, name);
    
    levV = linspace(rangeV(1),rangeV(2),levN);
    
    %mapV = makeColorMap([0 0 0.4; 0 0 0.8; 1 1 1; 1 0 0; 0.6 0 0], 200);
    %mapV = makeColorMap([0.4 0 0; 0.8 0 0; 1 1 1; 0 0 1; 0 0 0.6], 200);
    %mapV = makeColorMap([0 0 0; 0.15 0.15 1; 0.6 1 0.6; 0.9 1 0.9; 1 1 1],
    %200);
    %mapV = makeColorMap([0.1 0.1 0.2; 0.15 0.15 0.6; 0.15 1 0.6; 0.8 1 0.8; 1 1 1], 200);
    mapV = makeColorMap([0 0 0.4; 0.6 0.6 1; 1 1 1; 0.8 0.2 0.2; 1 0 0], 200);
    
    % Read the data.
    v = dlmread(crossFile, ' ');
    x = dlmread(crossFileX, ' ');
    y = dlmread(crossFileY, ' ');
    z = dlmread(crossFileZ, ' ');
      
    % Convert to nanometers from angstroms.
    z = 0.1*z;
    x = 0.1*x;
    y = 0.1*y;
    
    v = (v > 1e6).*1e6 + (v <= 1e6).*v;
    
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
    set(h, 'Position', [10 10 600 600], 'Resize', 'off');
    pcolor(cenX, cenY, cenV)
    %shading('interp')
    shading('flat')
    %contourf(cenX, cenY, cenV, levV)
    colormap(mapV)
    xlabel('x (nm)', 'FontSize', 20)
    ylabel('y (nm)', 'FontSize', 20)
    axis([-0.5*lFactor*lx 0.5*lFactor*lx -0.5*lFactor*ly 0.5*lFactor*ly])
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
    
    saveas(gcf, strcat(crossPlot,'.png'));
    print('-dpdf', strcat(crossPlot,'.pdf'));
end
