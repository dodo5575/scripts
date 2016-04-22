clear all;

binNum = 60;
%fileList = [dir('free*.dat'); dir('bound*.dat')]; 
fileList = [dir('free*timeRes2*raw*.dat')]; 

plotColor = [0 0 0; 0.9 0 0; 0.8 0.4 0; 0.5 0.5 0; 0 0.8 0; 0 0.7 0.7; 0 0 0.9; 0.6 0 0.6; 0.6 0.6 0.6];

figure(1)
clf
hold on

fprintf('\n****%s\n', datestr(now));
for j=1:length(fileList)
  fileName = fileList(j).name;
  data = dlmread(fileName, ' ');

  n = length(data(:,1));
  fprintf('\n%s: found %d items\n', fileName, n);

  % Data format: timeAvg xAvg dt dx
  tm = data(:,1);
  xm = data(:,2);
  dt = data(:,3);
  dx = data(:,4);

  fprintf('mean(dt): %g\n', mean(dt));
  fprintf('std(dt): %g\n', std(dt));
  fprintf('mean(dx): %g\n', mean(dx));
  fprintf('std(dx): %g\n', std(dx));
  
  % Histogram.
  [tauCount tauBin] = hist(dt, binNum);
  delTau = tauBin(2)-tauBin(1);
  tauProb = tauCount/(delTau*n);
  tauFile = sprintf('histTau_%s', fileName);
  dlmwrite(tauFile, [tauBin' tauProb'], ' ');
  
  [zetaCount zetaBin] = hist(dx, binNum);
  delZeta = zetaBin(2)-zetaBin(1);
  zetaProb = zetaCount/(delZeta*n);
  zetaFile = sprintf('histZeta_%s', fileName);
  dlmwrite(zetaFile, [zetaBin' zetaProb'], ' ');
  
  % Plot.
  color0 = plotColor(mod(j-1,length(plotColor))+1,:);
  gh = plot(dt, dx);
  set(gh, 'Color', color0);
  set(gh, 'Marker', 'o');
  set(gh, 'MarkerSize', 7);
  set(gh, 'LineStyle', 'None');  
end

hold off
