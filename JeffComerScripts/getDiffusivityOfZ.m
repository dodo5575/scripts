clear all;


% Paremeters:
segment0 = 11;
segmentWidth = 1;
cutTime = 0.05;
dimensions = 2;
% Input:
inDir = '.';
fileNum = 1:25;
filePre = 'dev_DoubleCut.0.';
fileSuf = '.msd';
% Output:
outFile = strcat('diffuse_', filePre,'dat');

plotColor = [0 0 0; 0.9 0 0; 0.8 0.4 0; 0.5 0.5 0; 0 0.8 0; 0 0.7 0.7; 0 0 0.9; 0.6 0 0.6; 0.6 0.6 0.6];
figure(1)
clf
hold on
fprintf('\n****%s\n', datestr(now));
out = fopen(outFile, 'w');

legendList={};
count = 0;
for j=1:length(fileNum)
    fileName = sprintf('%s/%s%d%s', inDir, filePre, fileNum(j), fileSuf);    
    % Get the center of the segment bin.
    segmentPos = segment0 + segmentWidth*fileNum(j);
    
    data = dlmread(fileName, ' ');
    if isempty(data)
       continue
    end
    dt = data(2,1) - data(1,1);
    cutIndex = ceil(cutTime/dt)+1;
    
    if length(data(:,1)) < cutIndex
        continue
    end
    
    % Add the legend entries.
    count = count + 1;
    legendList{2*(count-1)+1} = fileName;
    legendList{2*(count-1)+2} = fileName;
    
    % Get the data.
    t = data(cutIndex:end,1);
    n = length(t);
    msd = data(cutIndex:end,2);
    
    % Compute the linear regression.
    [a b da db] = linearRegression(t, msd);
    theo = b*t + a;
    
    % Compute the diffusivity.
    diffusion = b/(2*dimensions);
    diffusionErr = db/(2*dimensions);
    diffErr = 2*diffusionErr;
        
    % Plot it.
    color0 = plotColor(mod(j-1,length(plotColor))+1,:);
    gh = plot(t, msd);
    set(gh, 'Color', color0);
    set(gh, 'Marker', 'o');
    set(gh, 'MarkerSize', 7);
    set(gh, 'LineStyle', 'None');
   
    gh = plot(t, theo);
    set(gh, 'Color', color0);
    set(gh, 'LineStyle', '-');
    
    % Say the result.
    fprintf('diffusivity: %.3f %.3f +/- %.3f\n', segmentPos, diffusion, diffErr);
    fprintf(out, '%.10g %.10g %.10g\n', segmentPos, diffusion, diffErr);
end

fclose(out);

lh = legend(legendList);
set(lh, 'Interpreter', 'none');
xlabel('t (ns)');
ylabel('<x^2 + y^2> A^2');
hold off

