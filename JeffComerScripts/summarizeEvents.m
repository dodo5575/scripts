clear all;

globList = {'event_*dat'};

plotColor = [0 0 0; 0.9 0 0; 0.8 0.4 0; 0.5 0.5 0; 0 0.8 0; 0 0.7 0.7; 0 0 0.9; 0.6 0 0.6; 0.6 0.6 0.6];

fprintf('\n****%s\n', datestr(now));
for j=1:length(globList)
    fileList = dir(globList{j});
    fprintf('\nFound %d files for %s\n', length(fileList), globList{j});
    
    data = [];
    for f=1:length(fileList)
        d = dlmread(fileList(f).name, ' ');
        data = [data; d];
    end
  
    fprintf('samples %d\n', length(data(:,1)))
    fprintf('enterBoundFrac %g\n', mean(data(:,1)));
    fprintf('exitBoundFrac %g\n', mean(data(:,2)));
    fprintf('meanPassageTime %g\n', mean(data(:,3)));
    fprintf('meanBoundTime %g\n', mean(data(:,4)));
    fprintf('meanBoundDisp %g\n', mean(data(:,5)));
    fprintf('meanBoundNum %g\n', mean(data(:,6)));
    fprintf('meanFreeTime %g\n', mean(data(:,7)));
    fprintf('meanFreeDisp %g\n', mean(data(:,8)));
    fprintf('meanFreeNum %g\n', mean(data(:,9)));
        
    if (0)
        figure(1)
        clf
        plot(data(:,4), data(:,7), 'o');
        xlabel('exitTime (ns)');
        ylabel('boundMeanTime (ns)');
        
        figure(2)
        clf
        plot(data(:,4), data(:,9), 'o');
        xlabel('exitTime (ns)');
        ylabel('freeMeanTime (ns)');
    end
    
    % Histogram.
    %tau = data(:,4);
    %[histCount histTau] = hist(tau, 20);
    %dt = histTau(2)-histTau(1);
    %histProb = histCount/(n*dt);
    
    % Plot.
    %color0 = plotColor(mod(j-1,length(plotColor))+1,:)
    %gh = plot(histTau, histProb);
    %set(gh, 'Color', color0);
    %set(gh, 'Marker', 'none');
    %set(gh, 'MarkerSize', 7);
    %set(gh, 'LineStyle', '-');
end

