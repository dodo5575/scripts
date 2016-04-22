clear all;
dt = 1.0;
inFile = 'pore1.6_histo_data.txt';
outPrefix = 'histo1/';
%outPrefix = '';
dataDir = 'histo/';
xPrefix = 'cen_';
yPrefix = 'curr_';
%xLabel = 'COM position of res 127 (nm)';
xLabel = 'I/I_0';
yLabel = 'probability';
plotCurrent = 1;
ind = readIndex(inFile);
posMin = 2.0;
posMax = 4.0;

% Set open currents at various voltages.
current0 = [4.0 4.7 0.3; 6.5 7.6 0.5; 13 15.2 1.0];
plotColor = [0 0 0; 0.9 0 0; 0.8 0.4 0; 0.5 0.5 0; 0 0.8 0; 0 0.7 0.7; 0 0 0.9; 0.6 0 0.6; 0.6 0.6 0.6];

figure(1)
clf
hold all

outData = [];
n = length(ind(:,1));
for j=1:n
    name = ind{j,1};
    voltage = str2double(ind{j,2});
    offset = str2double(ind{j,4});
      
    currFactor = 1.0;
    for k=1:length(current0)
        if abs(voltage-current0(k,1)) < 0.25
            currFactor = 1.0/current0(k,2);
            break
        end
    end
    
    x = dlmread(sprintf('%s/%s%s.dat',dataDir,xPrefix,name), ' ');
    y = dlmread(sprintf('%s/%s%s.dat',dataDir,yPrefix,name), ' ');
    
    % Shift time by the offset.
    x(:,1) = x(:,1) + offset;
    y(:,1) = y(:,1) + offset;
    nt = length(y(:,1));
    
    count = 1;
    curr = [];
    for t=1:nt
        tNow = y(t,1);
        k = binsearch(tNow, x(:,1));
        
        % Skip if the time is not surrounded by position samples.
        if k<=0 || k>=length(x(:,1))
            continue
        end
        
        % Interpolate the position.
        t0 = x(k,1);
        t1 = x(k+1,1);
        pos0 = x(k,2);
        pos1 = x(k+1,2);
        pos = pos0 + (pos1-pos0)/(t1-t0)*(tNow-t0);

        % Select a section of the current samples.
        if pos >= posMin && pos < posMax
            curr(count) = y(t,2)*currFactor;
            count = count + 1;
        end
    end

    % Display the statistics.
    disp(sprintf('\n%s', name))
    disp(sprintf('current samples: %d', length(curr)))
    disp(sprintf('current mean: %f', mean(curr)))
    disp(sprintf('current std: %f', std(curr)))
    disp(sprintf('current se: %f', std(curr)/sqrt(length(curr))))
    
    % Write the statistics.
    fileName = strcat(outPrefix,'stat_',name,'.dat');
    fid = fopen(fileName,'wt');
    fprintf(fid, 'current samples: %d\n', length(curr));
    fprintf(fid, 'current mean: %f\n', mean(curr));
    fprintf(fid, 'current std: %f\n', std(curr));
    fprintf(fid, 'current se: %f\n', std(curr)/sqrt(length(curr)));
    fclose(fid);
        
    % Histogram.
    [num val] = hist(curr, 10);
    area = sum(num.*(val(2)-val(1)));
    prob = num/area;
        
    % Plot.
    if strcmp(ind{j,1}(1), 'l')
        shape = 's';
        color1 = 'k';
    elseif strcmp(ind{j,1}(1), 'c')
        shape = 'o';
        color1 = 'r';
    elseif strcmp(ind{j,1}(1), 'd')
        shape = 'd';
        color1 = 'b';
    end
    color0 = plotColor(mod(j-1,length(plotColor))+1,:);
    gh = plot(val,prob);
    set(gh, 'Color', color0);
    set(gh, 'Marker', shape);
    set(gh, 'MarkerSize', 7);
    set(gh, 'LineStyle', '-');
    
    % Write the results.
    if length(outPrefix) > 0
        dlmwrite(strcat(outPrefix,'histo_',yPrefix,name,'.dat'), [val' prob'], ' ');
    end
end

hold off
xlabel(xLabel)
ylabel(yLabel)
axis([-2 2 0 inf])
h = legend(ind{:,1});
set(h, 'Interpreter', 'none');



