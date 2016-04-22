clear all;

% Parameters:
outName = 'imgX_';
dir = 0;
run = 'cross/crossX_conc_POT_';
nameList = {'triplet_gcc' 'brown_triplet_gcc'};
        
rangeV = [0 6];
levN = 20;

lFactor = 1.0;
lx = 3.5;
ly = 6;
tickX = [-2 -1 0 1 2];
tickY = [-3 -2 -1 0 1 2 3];

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
    mapV = makeColorMap([0.15 0.15 0.2; 0.2 0.2 0.6; 0.2 0.7 0.6; ...
		    0.4 1 0.6; 0.7 1 0.8; 1 1 1], 200);
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
    set(h, 'Position', [10 10 300 600], 'Resize', 'off');
    pcolor(cenX, cenY, cenV)
    %shading('interp')
    shading('flat')
    %contourf(cenX, cenY, cenV, levV)
    colormap(mapV)
    xlabel('x (nm)', 'FontSize', 20)
    ylabel('y (nm)', 'FontSize', 20)
    axis([-0.5*lFactor*lx 0.5*lFactor*lx -0.5*lFactor*ly 0.5*lFactor*ly])
    axis equal
    caxis(rangeV)
    colorbar

    set(gcf, 'color', 'white'); % sets the color to white
    set(gca, 'TickDir', 'out', 'XTick', tickX, 'YTick', tickY, 'TickLength', [0.02 0.02]);
    set(gca, 'FontSize', 20);
    
    saveas(gcf, strcat(crossPlot,'.png'));
    print('-dpdf', strcat(crossPlot,'.pdf'));
    close(gcf);
end
