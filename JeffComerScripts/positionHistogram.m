% Histogram the position data.
clear all;
inGlob = 'run3_*timestep160*.dat';
sysLen = 32;
binSize = 0.15;
bulkZ = 22;
refFileInd = 1;
theoFile = '../theo/theo_energy_high.dat';

binEdges = (0:binSize:(sysLen-1))';
binCen = binEdges(1:end-1) + 0.5*binSize;

fprintf('\n****%s\n', datestr(now));
plotColor = [0 0 0; 0.9 0 0; 0.8 0.4 0; 0.5 0.5 0; 0 0.8 0; 0 0.7 0.7; 0 0 0.9; 0.6 0 0.6; 0.6 0.6 0.6];
figure(1)
clf
hold on

% Plot the files.
fileList = dir(strcat(inGlob));
nameList = cell(length(fileList), 1);
for f=1:length(fileList)
    nameList{f} = fileList(f).name;
    data = dlmread(fileList(f).name, ' ');
    n = length(data);
    
    count = histc(abs(data(:,3)), binEdges);
    count = count(1:end-1);
    den = count/n;
    
    if f == refFileInd
        refDen = den;
    end
    
    c = plotColor(mod(f-1,length(plotColor))+1,:);
    gh = plot(binCen, den);
    set(gh, 'Color', c);
    set(gh, 'Marker', 'o');
    set(gh, 'MarkerSize', 7);
    set(gh, 'LineStyle', '-');
end

theoData = dlmread(theoFile, ' ');
theoZ = theoData(:,1);
theoV = theoData(:,2);

% Find the
for j=1:length(binCen)
    if binCen(j) > bulkZ
        bulkInd = j;
        break;
    end
end
bulkConc = mean(refDen(bulkInd:end));
theoConc = bulkConc*exp(-theoV);
plot(theoZ, theoConc, 'b--');

hold off
xlabel('z (AA)');
ylabel('count');
axis([14 22 0 inf]);
h = legend(nameList);
set(h, 'Interpreter', 'none');
