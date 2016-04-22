clear all;
dt = 0.1;
inFile = 'index_curv.txt';
outFile = 'pos_vs_curv.dat';
% name(1) voltage(2) COM_position(3) current(4) position(5) stretch(6)
% angle(7)
xPick = 3;
yPick = 6;
xLabel = 'COM position (nm)';
yLabel = 'helix length (nm)';

ind = readIndex(inFile);

% Set open currents at various voltages.
current0 = [1.0 2.7 0.3; 1.5 3.9 0.1; 2.0 4.8 0.3; 2.5 5.0 1.0; ...
    3.0 5.305 1.0; 4.0 5.81 0.08; 6.0 7.5 0.1];
plotColor = [0 0 0; 0.9 0 0; 0.8 0.4 0; 0.5 0.5 0; 0 0.8 0; 0 0.7 0.7; 0 0 0.9; 0.6 0 0.6; 0.6 0.6 0.6];

figure(1)
clf
hold all

outData = [];
n = length(ind(:,1));
for j=1:n
    voltage = str2double(ind{j,2});
    
    currFactor = 1.0;
    for k=1:length(current0)
        if abs(voltage-current0(k,1)) < 0.25
            currFactor = 1.0/current0(k,2);
            break
        end
    end
    
    x = dlmread(ind{j,xPick}, ' ');
    y = dlmread(ind{j,yPick}, ' ');
       
    t = max([x(1,1) y(1,1)]);
    if x(end,1) < y(end,1)
        tf = x(end,1);
        nSamples = floor((tf-t)/dt);
    else
        tf = y(end,1);
        nSamples = floor((tf-t)/dt);
    end
       
    meanX = zeros(nSamples,1);
    seX = zeros(nSamples,1);
    meanY = zeros(nSamples,1);
    seY = zeros(nSamples,1);
    k = 1;
    while t < tf
        [mX sX nX] = intervalMean(t, t+dt, x(:,1), x(:,2));
        [mY sY nY] = intervalMean(t, t+dt, y(:,1), y(:,2));
        
        if nX > 2 && nY > 2
            meanX(k) = mX;
            seX(k) = sX;
            meanY(k) = mY;
            seY(k) = sY;
            k = k + 1;
        end
        t = t + dt;
    end
    
    if k < nSamples
        meanX = meanX(1:k);
        seX = seX(1:k);
        meanY = meanY(1:k);
        seY = seY(1:k);
    end
    
    if strcmp(ind{j,1}(1), 'l')
        shape = 's';
    elseif strcmp(ind{j,1}(1), 'c')
        shape = 'o';
    elseif strcmp(ind{j,1}(1), 'd')
        shape = '*';
    end
    c = plotColor(mod(j-1,length(plotColor))+1,:);
    %gh = errorbar(meanX, meanY, seY);
    gh = plot(meanX, meanY);
    set(gh, 'Color', c);
    set(gh, 'Marker', shape);
    set(gh, 'MarkerSize', 7);
    set(gh, 'LineStyle', '-');
    
    %outData = [outData; [meanX meanY seX seY]];
end

hold off
xlabel(xLabel)
ylabel(yLabel)
h = legend(ind{:,1});
set(h, 'Interpreter', 'none');

% Write the results.
dlmwrite(outFile, outData, ' ');



